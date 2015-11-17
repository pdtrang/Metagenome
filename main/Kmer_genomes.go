// example how to count k-mer frequencies in a set of reads, using multiple goroutines.

package main

import (
   "github.com/vtphan/kmers"
   "os"
   "bufio"
   "runtime"
   "sync"
   "math"
   "encoding/csv"
   "strconv"
   "io/ioutil"
   "fmt"
)


func CountFreq(readFile string, freq map[int]int, K int) {

   // Get all reads into a channel
   reads := make(chan string)
   go func() {
      f, err := os.Open(readFile)
      if err != nil {
         panic("Error opening " + readFile)
      }
      defer f.Close()
      scanner := bufio.NewScanner(f)
      i := 0.0
      for scanner.Scan() {
        if math.Mod(i, 2.0) == 1 {
         reads <- (scanner.Text())
        }
        i++
      }
      close(reads)
   }()

   numCores := runtime.NumCPU()
   //numCores := 1
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
            // Count2 counts on both strands. Count1 counts only on the main strand.
            c.Count2([]byte(read))
         }
      }()
   }

   // Finish counting
   wg.Wait()

   // Print out the result
/*   for kmer := range(freq) {
      fmt.Println(kmer,kmers.NumToKmer(kmer,K),freq[kmer])
   }*/

}


func error_check(err error) {
   if err != nil {
      panic(err)
   }
}
func getGSM(GSMfile string) ([]int){

   csvfile, err := os.Open(GSMfile)

   error_check(err)

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
      error_check(err)
      //fmt.Println(i)
      index = append(index, i)
    }

   return index
}

func init_freq(index []int, freq map[int]int){

   for k := 0; k < len(index); k++ {
         freq[index[k]] = 0
   }

}

func assignTocsv(Kmer_freq map[int][]int, index []int, freq map[int]int, result string, name map[string]int) {
   resultfile, err := os.Create(result+".csv")
   error_check(err)

   indexfile, err := os.Create(result+"_index.csv")
   error_check(err)

   rw_idx := csv.NewWriter(indexfile)
   idx := make([]string, 1)

   rw := csv.NewWriter(resultfile)
   head := make([]string, len(name)+1)
   head[0] = "kmer"

   for n := range name{
      head[name[n]+1] = n
   }
   
   returnError := rw.Write(head)
   error_check(returnError)
   rw.Flush()

   /*csvfile, err := os.Open(GSMfile)
   error_check(err)
   defer csvfile.Close()
   reader := csv.NewReader(csvfile)
   _, err = reader.Read()

   for {
      line, err := reader.Read()
      if err != nil {
          break
      } else {
         head[0] = line[0]
         kmer, _ := strconv.Atoi(line[0])
         // fmt.Println(kmers.NumToKmer(kmer,14))
         head[1] = strconv.Itoa(freq[kmer])
         returnError := rw.Write(head)
         error_check(returnError)
      }
   }*/
   
   for k:= range Kmer_freq {
      head[0] = strconv.Itoa(k)
      idx[0] = strconv.Itoa(k)
      
      for n := range name {
         head[name[n]+1] = strconv.Itoa(Kmer_freq[k][name[n]])
      }
      returnError := rw.Write(head)
      error_check(returnError)
      returnError = rw_idx.Write(idx)
      error_check(returnError)
   }

   /*for i := 0; i<len(index); i++ {
      head[0] = strconv.Itoa(index[i])
      head[1] = strconv.Itoa(freq[index[i]])
      returnError := rw.Write(head)
      error_check(returnError)
   }*/

   rw.Flush()
   rw_idx.Flush()
}

func addResult(freq map[int]int, idx int, Kmer_freq map[int][]int) {

   for k := range freq {
         Kmer_freq[k][idx] = freq[k]
   }

}

func main() {
   if len(os.Args) != 5 {
      panic("Must provide readfile and GSM from genomes, kmer length and read GSM result file name!")
   }

   files, _ := ioutil.ReadDir(os.Args[1])
   name := make(map[string]int)

   for i, fi := range files {
      name[fi.Name()] = i
   }

   //fmt.Println(name)
   //fmt.Println(len(name))

   var index []int
   index = getGSM(os.Args[2])
   //fmt.Println(len(index))
   
   K, _ := strconv.Atoi(os.Args[3])
   freq := make(map[int]int)

   Kmer_freq := make(map[int][]int)

   for k := 0; k < len(index); k++ {
      for i:= 0; i < len(name); i++ {
         Kmer_freq[index[k]] = append(Kmer_freq[index[k]], 0)
      }
   }

   //fmt.Println("len = ", len(Kmer_freq)) 

   fmt.Println("K = ", K)
   fmt.Println("Number of kmer to count = ", len(index))
   fmt.Println("Number of genomes = ", len(files))

   for _, fi := range files {
      fmt.Println("Processing : ", fi.Name())
      init_freq(index, freq)
      CountFreq(os.Args[1] + "/" + fi.Name(), freq, K)
      addResult(freq, name[fi.Name()], Kmer_freq)
   } 
   
   assignTocsv(Kmer_freq, index, freq, os.Args[4], name)
}
