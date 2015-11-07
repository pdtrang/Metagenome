/*
  Convert number to k-mer.
*/

package main

import (
   "fmt"
   "encoding/csv"
   "strconv"
   "os"
   "sync"
   "time"
   "io/ioutil"
   "runtime"
   "strings"
   "bufio"
)

/*
   Return the K-mer (consisting of A,C,G,T) represented by x
*/
func NumToKmer(x int, K int) string {
   y := make([]byte, K)
   for i:=K-1; i>=0; i-- {
      base := x%4
      switch base {
         case 0: y[i] = 'A'
         case 1: y[i] = 'C'
         case 2: y[i] = 'G'
         case 3: y[i] = 'T'
      }
      x = (x-base)>>2
   }
   return string(y)
}

func main(){
   files, _ := ioutil.ReadDir(os.Args[1])
   
   start_time := time.Now()

   var wg sync.WaitGroup
   namesChan := make(chan string, len(files))
   for idx, f := range(files) {
      namesChan <- (strconv.Itoa(idx) + "," + f.Name())
   }
   close(namesChan)

   workers := runtime.NumCPU()
   runtime.GOMAXPROCS(workers)
   fmt.Println()
   if len(files) < workers {
      workers = len(files)
   }

   file1, err := os.Create("makeRead")
   if err != nil {
      fmt.Println("error created file")
   }
   defer file1.Close()

   w1 := bufio.NewWriter(file1)

   //var mutex = &sync.Mutex{}
   for i := 0; i < workers; i++ {
      wg.Add(1)
      go func(){
         defer wg.Done()
         for name := range namesChan {
            elements := strings.Split(name, ",")
            //idx, _ := strconv.Atoi(elements[0])
            filename := "kmers_"+elements[1]
            filename = strings.Replace(filename, "fsa.csv", "fna", -1)
            
            fmt.Println(filename)
            filefq := strings.Replace(filename, ".fna", ".fq", -1)
            fmt.Fprint(w1, "python2 make_real_read.py output_test/",filename, " output_test/realReads/",filefq,"\n")

            csvfile, err := os.Open(os.Args[1]+elements[1])

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

            K:= 14
            var kmer []string
            for i:=0; i< len(index);i++ {
               y:= NumToKmer(index[i],K)
               kmer = append(kmer, y)
            }
                                    
            file, err := os.Create(os.Args[2]+filename)
            if err != nil {
               fmt.Println("error created file")
            }
            defer file.Close()

            w := bufio.NewWriter(file)

            for i:=0; i<len(kmer);i++ {
               fmt.Fprint(w, ">k"+strconv.Itoa(index[i]), "\n")
               fmt.Fprint(w, kmer[i], "\n")
            }
           
            w.Flush()

            w1.Flush()
       
            
         }
      }()
   }
   wg.Wait()
   
   end_time := time.Since(start_time)
   fmt.Println("Used time", end_time)
}
