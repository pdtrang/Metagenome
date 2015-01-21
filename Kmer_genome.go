package main

import (
   "os"
   "fmt"
   "bufio"
   "bytes"
   "sync"
   "runtime"
   "strconv"
   "time"
   "io/ioutil"
   "encoding/csv"
   //"sort"
)

// os.Argv[1]: sequences folder
// os.Argv[2]: result filename

func main() {
    if len(os.Args) != 3 {
        panic("must provide sequence folder file and result file name.")
    }

    files, _ := ioutil.ReadDir(os.Args[1])
    start_time := time.Now()
    gsm := make(map[int][]int)
    gsmFreq := make(map[int][]int)

    core_num := 2
    kmer_len := 7
    distance := 0
    runtime.GOMAXPROCS(core_num+2)

    resultfile, err := os.Create(os.Args[2]+".csv")
    if err != nil {
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }
    rw := csv.NewWriter(resultfile)
    head := make([]string, len(files)+1)
    head[0] = "K-mer"

    for index, fi := range files {
        head[index+1] = fi.Name()
    }

    returnError := rw.Write(head)
    if returnError != nil { 
        fmt.Println(returnError)
    }
    rw.Flush()

    for index, fi := range files {
        f,err := os.Open(os.Args[1] + "/" + fi.Name())
        if err != nil {
            fmt.Printf("%v\n",err)
            os.Exit(1)
        }

        br := bufio.NewReader(f)
        byte_array := bytes.Buffer{}

        _, isPrefix, err := br.ReadLine()

        if err != nil || isPrefix {
            fmt.Printf("%v\n",err)
            os.Exit(1)
        }

        for {
            line , isPrefix, err := br.ReadLine()
            if err != nil || isPrefix{
                break
            } else {
                byte_array.Write([]byte(line))
            }       
        }

        input := []byte(byte_array.String())
        var wg sync.WaitGroup
        result := make(chan int, core_num)

        for i := 0; i < core_num; i++ {
            wg.Add(1)
            go process(input, i, core_num, kmer_len, distance, result, &wg)
        }

        go func() {
            wg.Wait()
            close(result)
        }()

        gsm1 := make(map[int]int)
        gsmFreq1 := make(map[int]int)
        
        fmt.Println("processing : ", fi.Name())
        
        for res := range result {
            gsm1[res] = index + 1
            gsmFreq1[res] = gsmFreq1[res] + 1
        }

        for k := range gsm1 {
                gsm[k] = append(gsm[k], gsm1[k])
                gsmFreq[k] = append(gsmFreq[k], gsmFreq1[k])
        }
        
        f.Close()

    }
    
    //save index
    idx, err := os.Create("index.csv")
    if err != nil {
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }
    rrk := csv.NewWriter(idx)

    for k := range gsm {
        
        line := make([]string, len(files)+1)
        line_idx := make([]string, 1)
        for i := range line {
            if i == 0 {
                line[0] = strconv.Itoa(k)
                line_idx[0] = strconv.Itoa(k)
            } else {
                line[i] = strconv.Itoa(0)
                for j:=0; j<len(gsm[k]); j++ {
                    if i == gsm[k][j] {
                        line[i] = strconv.Itoa(gsmFreq[k][j])
                    } 
                }
            }   

        }
        
        returnError1 := rw.Write(line)
        if returnError1 != nil {
            fmt.Println(returnError1)
        }
        returnError_idx := rrk.Write(line_idx)
        if returnError_idx != nil {
            fmt.Println(returnError_idx)
        } 
        
    }
    rw.Flush()
    rrk.Flush()

    gsm_time := time.Since(start_time)
    fmt.Println("used time", gsm_time)
}

func process(genome []byte, i int, core_num int, kmer_len int, distance int, result chan int, wg *sync.WaitGroup) {
    defer wg.Done()
    begin := len(genome)*i/core_num
    end := len(genome)*(i+1)/core_num
    if begin != 0 {
        begin = begin - kmer_len + 1
    }
    for m := begin; m < end-kmer_len+1; m++ {
        m1 := m
        m2 := m+kmer_len
        kmer := make([]byte, kmer_len)
        copy(kmer, genome[m1:m2])

        if (distance > 0){
            m3 := m+kmer_len+distance
            m4 := m+2*kmer_len+distance    
            kmer = append(kmer, genome[m3:m4]...)
       }
               
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
        
        if repr!= -1 {
            result <- repr
        }
    }
    
}

