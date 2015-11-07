package main

import (
   "os"
   "fmt"
   "time"
   "strconv"
   "encoding/csv"
   "io/ioutil"
   "strings"
   "math"
   "sync"
   "runtime"
//   "bufio"
)

//parallel
//os.Args[1]: folder for F
//os.Args[2]: number of gsm
//os.Args[3]: proportion
//os.Args[4]: output file

func Read_total_gsm() (map[string]int, map[string]int){
   csvfile, err := os.Open(os.Args[2])

   if err != nil {
      fmt.Println(err)
      return nil, nil
   }

   defer csvfile.Close()
   reader := csv.NewReader(csvfile)
   reader.FieldsPerRecord = -1 
   rawCSVdata, err := reader.ReadAll()

   if err != nil {
      fmt.Println(err)
      os.Exit(1)
   }
    
   //var name []string
   gsm := make(map[string]int)
   length := make(map[string]int)

   for _, each := range rawCSVdata {
      n := each[0]
      g, err := strconv.Atoi(each[1])
      l, err := strconv.Atoi(each[2])
      if err != nil {
         // handle error
         fmt.Println(err)
         os.Exit(2)
      }
      //fmt.Println(i)
      //name = append(name, n)
      gsm[n] = g
      length[n] = l
   }
   
   return gsm, length
}

func Read_Fc(readFile string) ([]int, []int){
   csvfile, err := os.Open(readFile)

   if err != nil {
      fmt.Println(err)
      return nil, nil
   }

   defer csvfile.Close()
   reader := csv.NewReader(csvfile)
   reader.FieldsPerRecord = -1 
   rawCSVdata, err := reader.ReadAll()

   if err != nil {
      fmt.Println(err)
      os.Exit(1)
   }
    
   var F []int
   var c []int
   
   for _, each := range rawCSVdata {
      f, err := strconv.Atoi(strings.Replace(each[1], " ", "", -1))
      b, err := strconv.Atoi(strings.Replace(each[2], " ", "", -1))
      if err != nil {
         // handle error
         fmt.Println(err)
         os.Exit(2)
      }
      //fmt.Println(i)
      F = append(F, f)
      c = append(c, b)   
            
   }
   
   return F, c
}

func Sort_count(F []int, c []int) {
   for i:=0; i<len(c)-1; i++ {
      for j:=i+1;j<len(c);j++ {
         if (c[i] < c[j]){
            temp := c[i]
            c[i] = c[j]
            c[j] = temp

            temp = F[i]
            F[i] = F[j]
            F[j] = temp
         }
      }
   }

}

// https://gist.github.com/DavidVaini/10308388
func Round(val float64, roundOn float64, places int ) (newVal float64) {
   var round float64
   pow := math.Pow(10, float64(places))
   digit := pow * val
   _, div := math.Modf(digit)
   if div >= roundOn {
      round = math.Ceil(digit)
   } else {
      round = math.Floor(digit)
   }
   newVal = round / pow
   return
}

// conpute n by genome length
func Compute_n_glength(kmer int, pr float64, gsm map[string]int, length map[string]int) (map[string]int){
   n := make(map[string]int)
   
   for k := range gsm {
      temp := int(Round((pr*float64(length[k]))/float64(kmer), 0.5, 0))
      //temp := int(Round((pr*gsm[k]), 0.5, 0))
      if (temp == 0){
         n[k] = 1
      }else{
         n[k] = int(math.Min(float64(temp), float64(gsm[k])))/2
      }
      //fmt.Println(k, n[k])
   }
   

   return n
}

// compute n by # of gsm
func Compute_n_gsm(kmer int, pr float64, gsm map[string]int, length map[string]int) (map[string]int){
   n := make(map[string]int)
   
   for k := range gsm {
      temp := int(Round((pr*float64(gsm[k])), 0.5, 0))
      if (temp == 0){
         n[k] = 1
      }else{
         n[k] = int(math.Min(float64(temp), float64(gsm[k])))/2
      }
      //fmt.Println(k, n[k])
   }
   

   return n
}

//If n is odd then Median (M) = value of ((n + 1)/2)th item term.
//If n is even then Median (M) = value of [((n)/2)th item term + ((n)/2 + 1)th item term ]/2
func Get_range_median(c []int, l int) (int, int){
   s := 0
   e := 0
   if (len(c) > 2 * l){
      if(len(c) % 2 == 1){
         s = (len(c)/2) - l
         e = (len(c)/2) + l + 1
      }else{
         s = (len(c)/2) -l
         e = (len(c)/2) + l 
      }
   }else{
      s = 0
      e = len(c)
   }

   return s, e
}

func SaveToVector(res map[string]float64, outfile string){
   //create output file
   output := outfile
   fmt.Println("Saving output: ", outfile)

   resultread, err := os.Create(output)
   if err != nil {
      fmt.Printf("%v\n",err)
      os.Exit(1)
   }
   rr := csv.NewWriter(resultread)
    
   //save to file
   for k := range res {
      line:= make([]string, 2)
      line[0] = k
      line[1] = strconv.FormatFloat(res[k], 'f', -1, 64)
      returnError := rr.Write(line)
      if returnError != nil {
         fmt.Println(returnError)
      }
   }
   rr.Flush()
}

/////////////////////////
// compute by n_length //
/////////////////////////
func MedianByLength(F []int, c []int, n_length int) (float64){
   s, e := Get_range_median(c, n_length)
         
   sum_F := 0
   sum_c := 0
   for i := s; i<e; i++ {
      sum_F = sum_F + F[i]
      sum_c = sum_c + c[i]
   }

   //fmt.Println("Median GSM by Prop to Genome Length: ", filename, n_length, s, e, sum_F, sum_c)
   
   r1 := float64(sum_c)/float64(sum_F)

   return r1
}

//////////////////////
// compute by n_gsm //
//////////////////////
func MedianByGSM(F []int, c []int, n_gsm int) (float64){
   s, e := Get_range_median(c, n_gsm)
         
   sum_F := 0
   sum_c := 0
   for i := s; i<e; i++ {
      sum_F = sum_F + F[i]
      sum_c = sum_c + c[i]
   }

   //fmt.Println("Median GSM by Prop to Total Genome GSM: ", filename, n_gsm, s, e, sum_F, sum_c)
   
   r2 := float64(sum_c)/float64(sum_F) 

   return r2
}

/////////////////////////////
// compute by top n_length //
/////////////////////////////
func TopByLength(F []int, c []int, n_length int) (float64){
   s, e := 0, 2*n_length
         
   sum_F := 0
   sum_c := 0
   for i := s; i<e; i++ {
      sum_F = sum_F + F[i]
      sum_c = sum_c + c[i]
   }

   //fmt.Println("Top GSM by Prop to Genome Length: ", filename, n_length, s, e, sum_F, sum_c)
   
   r3 := float64(sum_c)/float64(sum_F) 

   return r3 
}

//////////////////////////
// compute by top n_gsm //
//////////////////////////
func TopByGSM(F []int, c []int, n_gsm int) (float64){
   s, e := 0, 2*n_gsm
         
   sum_F := 0
   sum_c := 0
   for i := s; i<e; i++ {
      sum_F = sum_F + F[i]
      sum_c = sum_c + c[i]
   }

   //fmt.Println("Top GSM by Prop to Total Genome GSM: ", filename, n_gsm, s, e, sum_F, sum_c)
   
   r4 := float64(sum_c)/float64(sum_F)

   return r4
}

func Process(file string, n_length int, n_gsm int, filename string, idx int) (float64, float64, float64, float64){
   var F []int
   var c []int
   F, c = Read_Fc(file)    
   
   Sort_count(F, c)

   r1 := MedianByLength(F, c, n_length)
   r2 := MedianByGSM(F, c, n_gsm)
   r3 := TopByLength(F, c, n_length)
   r4 := TopByGSM(F, c, n_gsm)
   
   return r1, r2, r3, r4
}

func main() {
   if len(os.Args) != 5 {
        panic("must provide folder for F, total gsm file, proportion and output file name.")
   }

   files, _ := ioutil.ReadDir(os.Args[1])
   k := 14
   start_time := time.Now()
   
   gsm := make(map[string]int)
   length := make(map[string]int)
   n_length := make(map[string]int)
   n_gsm := make(map[string]int)
   gsm, length = Read_total_gsm()
   pr, err := strconv.ParseFloat(os.Args[3], 64)
   if err != nil {
      panic("pr error")
   }

   if (pr < 1.0){
      n_length = Compute_n_glength(k, pr, gsm, length)
      n_gsm = Compute_n_gsm(k, pr, gsm, length)
   }else{
      n_length = gsm
      n_gsm = gsm
   }

   res_median_length := make(map[string]float64)
   res_median_gsm := make(map[string]float64)
   res_top_length := make(map[string]float64)
   res_top_gsm := make(map[string]float64)
   res_length := make(map[string]float64)
   res_gsm := make(map[string]float64)
   
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
   
   var mutex = &sync.Mutex{}
   for i := 0; i < workers; i++ {
      wg.Add(1)
      go func(){
         defer wg.Done()
         for name := range namesChan {
            elements := strings.Split(name, ",")
            idx, _ := strconv.Atoi(elements[0])
            filename := elements[1]
            
            fmt.Println(filename)
            r1, r2, r3, r4 := Process(os.Args[1] + "/" + filename, n_length[filename], n_gsm[filename], filename, idx)

            mutex.Lock()
            res_median_length[filename] = r1
            res_median_gsm[filename] = r2
            res_top_length[filename] = r3
            res_top_gsm[filename] = r4
            res_length[filename] = (r1+r3)/2
            res_gsm[filename] = (r2+r4)/2
            mutex.Unlock()
            //res[filename] = float64(sum_c)/float64(sum_F)                     
         }
      }()
   }
   wg.Wait()

   output1 := "MedianByLength_"+os.Args[4]
   output2 := "MedianByGSM_"+os.Args[4]
   output3 := "TopByLength_"+os.Args[4]
   output4 := "TopByGSM_"+os.Args[4]
   output5 := "Length_"+os.Args[4]
   output6 := "GSM_"+os.Args[4]
   SaveToVector(res_median_length, output1)
   SaveToVector(res_median_gsm, output2)
   SaveToVector(res_top_length, output3)
   SaveToVector(res_top_gsm, output4)
   SaveToVector(res_length, output5)
   SaveToVector(res_gsm, output6)

   end_time := time.Since(start_time)
   fmt.Println("Used time", end_time)
}
