/*
Author: Vinhthuy Phan, 2015
*/
package Metagenome

import (
   "sync"
   "math"
   // "fmt"
)

/*
   Counter is used to counter ALL kmers in sequences and their reverse complements.
*/
type Counter struct {
   K int
   Freq []int
   lock []sync.RWMutex
}

func NewCounter(K int) *Counter {
   c := new(Counter)
   c.K = K
   c.Freq = make([]int, int(math.Pow(4,float64(K))))
   c.lock = make([]sync.RWMutex, len(c.Freq))
   return c
}

func (c *Counter) Count(sequence []byte) {
   K := c.K
   uK := uint(K)
   var id int
   var id_rc int
   for i:=0; i<len(sequence)-K+1; i++ {
      if i==0 || (sequence[i+K-1]!='A' && sequence[i+K-1]!='C' && sequence[i+K-1]!='G' && sequence[i+K-1]!='T') {
         i = firstIndex(sequence, K, i, len(sequence))
         if i==-1 {
            return
         }
         id, id_rc = Kmer2Num(sequence, c.K, i)
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
      c.lock[id].Lock()
      c.Freq[id]++
      c.lock[id].Unlock()

      c.lock[id_rc].Lock()
      c.Freq[id_rc]++
      c.lock[id_rc].Unlock()
   }
}

// ----------------------------------------------------------------

/*
   KmerCounter is used to counter only kmers in Freq in reads.  kmers in
   reverse complements of reads are *not* counted.
*/
type KmerCounter struct {
   K int
   Freq map[int]int
   lock map[int]*sync.RWMutex
}

/*
   Count all kmers in freq
*/
func NewKmerCounter(K int, freq map[int]int) *KmerCounter {
   c := new(KmerCounter)
   c.K = K
   c.Freq = freq
   c.lock = make(map[int]*sync.RWMutex)

   for kmer := range(freq) {
      c.lock[kmer] = new(sync.RWMutex)
   }
   return c
}

func (c *KmerCounter) Count(sequence []byte) {
   K := c.K
   uK := uint(K)
   var id int
   var ok bool
   for i:=0; i<len(sequence)-K+1; i++ {
      if i==0 || (sequence[i+K-1]!='A' && sequence[i+K-1]!='C' && sequence[i+K-1]!='G' && sequence[i+K-1]!='T') {
         i = firstIndex(sequence, K, i, len(sequence))
         if i==-1 {
            return
         }
         id, _ = Kmer2Num(sequence, c.K, i)
      } else {
         switch sequence[i-1] {
            case 'A':
               id = (id<<2)
            case 'C':
               id = (id<<2) - (1 << (uK<<1))
            case 'G':
               id = (id<<2) - (2 << (uK<<1))
            case 'T':
               id = (id<<2) - (3 << (uK<<1))
            default:
               panic("Count: invalid base " + string(sequence[i-1]))
         }
         switch sequence[i+K-1] {
            case 'A':
            case 'C':
               id += 1
            case 'G':
               id += 2
            case 'T':
               id += 3
            default:
               panic("Count: invalid base " + string(sequence[i+K-1]))
         }
      }
      _, ok = c.Freq[id]
      if ok {
         c.lock[id].Lock()
         c.Freq[id]++
         c.lock[id].Unlock()
      }
   }
}

