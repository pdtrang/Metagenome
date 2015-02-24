/*
Author: Vinhthuy Phan, 2014
Get all kmers of a sequence
*/
package Metagenome

import (
   "sync"
   // "fmt"
)

/*
   compute the (base-4) numeric representations of the sequence
   and its reverse complement
*/
func Kmer2Num(sequence []byte, K int, i int) (int, int) {
   id := 0
   id_rc := 0
   for j:=i; j<K+i; j++ {
      switch sequence[j] {
         case 'A': id = id<<2
         case 'C': id = id<<2 + 1
         case 'G': id = id<<2 + 2
         case 'T': id = id<<2 + 3
         default:
            return -1, -1
      }
      switch sequence[K+i-1+i-j] {
         case 'T': id_rc = id_rc<<2
         case 'G': id_rc = id_rc<<2 + 1
         case 'C': id_rc = id_rc<<2 + 2
         case 'A': id_rc = id_rc<<2 + 3
         default:
            return -1, -1
      }
   }
   return id, id_rc
}

/*
   The first index, starting from "start" forwardly, of a kmer consisting of
   only A,C,T,G.
*/
func firstIndex(sequence []byte, K int, start int, end int) int {
   var j int
   for i:=start; i<end-K+1; i++ {
      for j=i; j<i+K; j++ {
         if sequence[j]!='A' && sequence[j]!='C' && sequence[j]!='G' && sequence[j]!='T' {
            break
         }
      }
      if j==i+K {
         return i
      }
   }
   return -1
}


/*
Store in channel "result" all kmers (represented by numbers) of "sequence" and its
reverse complement by sliding a window of length K from "start" to "end"-1.

h := map[byte]int{'A':0, 'C':1, 'G':2, 'T':3}
id = (id<<2) - (h[sequence[i-1]] << (uint(2*K))) + h[sequence[i+K-1]]
*/
func Slide(sequence []byte, K int, start int, end int, result chan int) {
   uK := uint(K)
   var id int
   var id_rc int
   for i:=start; i<end-K+1; i++ {
      if i==start || (sequence[i+K-1]!='A' && sequence[i+K-1]!='C' && sequence[i+K-1]!='G' && sequence[i+K-1]!='T') {
         i = firstIndex(sequence, K, i, end)
         if i==-1 {
            return
         }
         id, id_rc = Kmer2Num(sequence, K, i)
      } else {
         switch sequence[i-1] {
            case 'A':
               id = (id<<2)
               id_rc = (id_rc - 3) >> 2
            case 'C':
               id = (id<<2) - (1 << (uK<<1))
               id_rc = (id_rc - 2) >> 2
            case 'G':
               id = (id<<2) - (2 << (uK<<1))
               id_rc = (id_rc - 1) >> 2
            case 'T':
               id = (id<<2) - (3 << (uK<<1))
               id_rc = id_rc >> 2
            default:
               panic("Slide: invalid base " + string(sequence[i-1]))
         }
         switch sequence[i+K-1] {
            case 'A':
               id_rc += 3 << ((uK-1)<<1)
            case 'C':
               id += 1
               id_rc += 2 << ((uK-1)<<1)
            case 'G':
               id += 2
               id_rc += 1 << ((uK-1)<<1)
            case 'T':
               id += 3
            default:
               panic("Slide: invalid base " + string(sequence[i-1]))
         }
      }
      result <- id
      result <- id_rc
   }
}

/*
   Count all k-mers in sequence and its reverse complement by sliding a window across the sequence
*/
func Count(sequence []byte, K int, start int, end int, freq []int, lock []sync.RWMutex) {
   uK := uint(K)
   var id int
   var id_rc int
   for i:=start; i<end-K+1; i++ {
      if i==start || (sequence[i+K-1]!='A' && sequence[i+K-1]!='C' && sequence[i+K-1]!='G' && sequence[i+K-1]!='T') {
         i = firstIndex(sequence, K, i, end)
         if i==-1 {
            return
         }
         id, id_rc = Kmer2Num(sequence, K, i)
      } else {
         switch sequence[i-1] {
            case 'A':
               id = (id<<2)
               id_rc = (id_rc - 3) >> 2
            case 'C':
               id = (id<<2) - (1 << (uK<<1))
               id_rc = (id_rc - 2) >> 2
            case 'G':
               id = (id<<2) - (2 << (uK<<1))
               id_rc = (id_rc - 1) >> 2
            case 'T':
               id = (id<<2) - (3 << (uK<<1))
               id_rc = id_rc >> 2
            default:
               panic("Count: invalid base " + string(sequence[i-1]))
         }
         switch sequence[i+K-1] {
            case 'A':
               id_rc += 3 << ((uK-1)<<1)
            case 'C':
               id += 1
               id_rc += 2 << ((uK-1)<<1)
            case 'G':
               id += 2
               id_rc += 1 << ((uK-1)<<1)
            case 'T':
               id += 3
            default:
               panic("Count: invalid base " + string(sequence[i+K-1]))
         }
      }
      // fmt.Println(i, id, "\n", NumToKmer(id, K), "\n", string(sequence[i:i+K]))
      lock[id].Lock()
      freq[id]++
      lock[id].Unlock()

      lock[id_rc].Lock()
      freq[id_rc]++
      lock[id_rc].Unlock()
   }
}


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
