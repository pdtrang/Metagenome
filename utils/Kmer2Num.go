package main

import (
  "fmt"
  "encoding/csv"
  "os"
  "strconv"
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

func Kmer2Num(sequence []byte, K int, i int) (int) {
    id := 0
    
    for j:=i; j<K+i; j++ {
        switch sequence[j] {
            case 'A': id = id<<2
            case 'C': id = id<<2 + 1
            case 'G': id = id<<2 + 2
            case 'T': id = id<<2 + 3
            default:
                return -1
        }
            }

   return id
}

func main(){
  //start_time := time.Now()

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

  var kmer []string
  
  for _, each := range rawCSVdata {
    i := each[0]
    
    //fmt.Println(i)
    kmer = append(kmer, i)
  }
  //finish reading index file

  K:= 10
  var index []int
  for i:=0; i< len(kmer);i++ {
    y:= Kmer2Num([]byte(kmer[i]),K, 0)
    index = append(index, y)
  }

  //fmt.Println(len(index))
  for i:=0; i< len(index);i++ {
    fmt.Println(strconv.Itoa(index[i])+",")
  }


  //end_time := time.Since(start_time)
  //fmt.Println("Used time", end_time)
}
