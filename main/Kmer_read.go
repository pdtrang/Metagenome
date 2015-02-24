package main

import (
   "github.com/vtphan/kmers"
   "os"
   "bufio"
   "fmt"
   "runtime"
   "sync"
   "time"
   "strconv"
   "strings"
   "encoding/csv"
)

//os.Args[1]: index file
//os.Args[2]: reads
//os.Args[3]: output folder

func CountFreq(readFile string, K int, freq map[int]int) {

   fmt.Println("\nProcessing: ", readFile)
   // Get all reads into a channel
   reads := make(chan string)
   go func() {
      f, err := os.Open(readFile)
      if err != nil {
         panic("Error opening " + readFile)
      }
      defer f.Close()
      scanner := bufio.NewScanner(f)
      for scanner.Scan() {
         reads <- scanner.Text()
      }
      close(reads)
   }()

   numCores := runtime.NumCPU()
   runtime.GOMAXPROCS(numCores)
   var wg sync.WaitGroup

   // Start a new counter that counts only kmers in freq.
   c := kmers.NewKmerCounter(K, freq)

   // Count kmers in different cores simultaneously.
   for i:=0; i<numCores; i++ {
      wg.Add(1)
      go func() {
         defer wg.Done()
         for read := range(reads){
            c.Count([]byte(read))
         }
      }()
   }

   // Finish counting
   wg.Wait()
}

func Read_index() ([]int){
   csvfile, err := os.Open(os.Args[1])

   if err != nil {
      fmt.Println(err)
      return nil
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
   
   return index
}

func SaveToVector(index []int, freq map[int]int){
   //create output file
   out := strings.Split(os.Args[2], "/")
   out[len(out)-1] = "read_" + out[len(out)-1] + ".csv"
   output := os.Args[3] + out[len(out)-1]
   fmt.Println("Vector b: ", out[len(out)-1])

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

   //save to file
   for k:= range (index){
      line:=make([]string, 2)
      line[0] = strconv.Itoa(index[k])
      line[1] = strconv.Itoa(freq[index[k]])
      returnError := rr.Write(line)
      if returnError != nil {
         fmt.Println(returnError)
      }
   }
   rr.Flush()
}

func main() {
   if len(os.Args) != 4 {
        panic("must provide index file, readfile and result file name.")
   }

   start_time := time.Now()
   
   var index []int
   index = Read_index()
         
   K := 7
   freq := make(map[int]int)
   for i := range(index){
      freq[index[i]] = 0
   }
   CountFreq(os.Args[2], K, freq)
 
   SaveToVector(index, freq)
   
   end_time := time.Since(start_time)
   fmt.Println("Used time", end_time)
}
