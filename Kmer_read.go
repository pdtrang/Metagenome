package main

import (
   "os"
   "fmt"
   "bufio"
   
   "sync"
   "runtime"
   "strconv"
   "time"
   "strings"
   "encoding/csv"
   
)

//os.Args[1]: index file
//os.Args[2]: reads
//os.Args[3]: output folder

func main() {
    if len(os.Args) != 4 {
        panic("must provide index file, readfile and result file name.")
    }

    start_time := time.Now()

    out := strings.Split(os.Args[2], "/")
    resultRead := make(chan int, 10)

    //count kmer in reads
    go CountFreq(os.Args[2], 7, resultRead)
    kmrread := make(map[int]int)
    for res := range resultRead {
        kmrread[res] = kmrread[res] + 1
    }
    //fmt.Println(kmrread)

    //begin reading index file
    csvfile, err := os.Open(os.Args[1])

    if err != nil {
         fmt.Println(err)
         return
    }

    defer csvfile.Close()

    reader := csv.NewReader(csvfile)

    reader.FieldsPerRecord = -1 

    rawCSVdata, err := reader.ReadAll()

    if err != nil {
         fmt.Println(err)
         os.Exit(1)
    }

    
    var index []int
    for _, each := range rawCSVdata {
        i, err := strconv.Atoi(each[0])
        if err != nil {
            // handle error
            fmt.Println(err)
            os.Exit(2)
        }
        //fmt.Println(i)
        index = append(index, i)
    }
    //finish reading index file
 
    //create output file
    out[len(out)-1] = "read_" + out[len(out)-1] + ".csv"
    output := os.Args[3] + out[len(out)-1]
    fmt.Println(out[len(out)-1])
    resultread, err := os.Create(output)
    if err != nil {
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }
    rr := csv.NewWriter(resultread)
    head1 := make([]string, 2)
    head1[0] = "K-mer"
    head1[1] = "b"
    
    returnError1 := rr.Write(head1)
    if returnError1 != nil { 
        fmt.Println(returnError1)
    }
    rr.Flush()

    //save kmrread to output file
    for k := range index {
        key := index[k]
        if kmrread[k] != -1 {
            line := make([]string, 2)
            for i := range line {
                if i == 0 {
                    line[0] = strconv.Itoa(key)
                } else if i == 1 {
                    line[i] = strconv.Itoa(kmrread[key])
                } else {
                    line[i] = strconv.Itoa(0)
                } 
            }
            returnError := rr.Write(line)
            if returnError != nil {
                fmt.Println(returnError)
            }
        }
    }
    rr.Flush()

    gsm_time := time.Since(start_time)
    fmt.Println("used time", gsm_time)
}

func CountFreq(readFile string, K int, result chan int) {

   // Get all reads into a channel
   reads := make(chan []byte)
   go func() {
      f, err := os.Open(readFile)
      if err != nil {
         panic("Error opening " + readFile)
      }
      defer f.Close()
      scanner := bufio.NewScanner(f)
      
      for scanner.Scan() {
        
        reads <- []byte(scanner.Text())   
        
      }
      close(reads)
   }()

   // Spread the processing of reads into different cores
   numCores := runtime.NumCPU()
   runtime.GOMAXPROCS(numCores)
   var wg sync.WaitGroup

   for i:=0; i<numCores; i++ {
      wg.Add(1)
      go func() {
         defer wg.Done()
         ProcessRead(reads, K, result)
      }()
   }

   go func() {
        wg.Wait()
        close(result)
   }()
}


func ProcessRead(reads chan []byte, kmer_len int, result chan int) {
    for read := range reads {
        for m := 0; m < len(read) - kmer_len +1; m++ {
            m1 := m
            m2 := m+kmer_len

            kmer := make([]byte, kmer_len)
            copy(kmer, read[m1:m2])
                      
            //fmt.Println(kmer)
            repr := 0
            d:
            for j := 0; j<len(kmer); j++ {
                switch kmer[j] {
                    case 'A': repr = 4*repr
                    case 'C': repr = 4*repr + 1
                    case 'G': repr = 4*repr + 2
                    case 'T': repr = 4*repr + 3
                    default:
                    // we skip any qgram that contains a non-standard base, e.g. N
                      repr = -1
                      break d
                }

            }
            if repr != -1 {
              result <- repr
            }
        }
    }
}
