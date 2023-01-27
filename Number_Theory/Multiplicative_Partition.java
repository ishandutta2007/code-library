package com.euler;

import java.io.PrintStream;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import com.euler.common.EulerUtils;
import com.euler.common.EulerUtils.Pair;
import com.euler.common.Primes;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.SetMultimap;

public class Euler452 {
    private final static int N=1000000000;
    private final static long MOD=1234567891l;
    /*
     * All over this program, we will identify abstract products using prime indices repetitions.
     *  * The class "factor decomposition" is just a (descending) ordered list of prime factors.
     *      For example, for 12 we have [2,1]. For 30 we have [1,1,1]. For 16 we have [4]. For
     *      every prime number, it's just [1].
     *      This class is used later as a holder of non-necessarily prime factors.
     *     * The class "abstract product" holds a map, from prime index to amount of times.
     *         That is, if the abstract products is the divisor of an abstract number,
     *         N = p1^3*p2^2*p3^2, the product may have 0, 1 or 2 as prime indices (map key),
     *         and the value can't be (respectively) higher than 3, 2 or 2.
     *  * The class "multiplicative partition" holds an indefinite amount of non-empty
     *      abstract products, so that all of them amount to the original number.
     *      Following the previous example, several possible partitions are:
     *          {[0->3, 1->2, 2->2]}
     *          {[0->1], [0->1], [0->1, 1->2, 2->2]}
     *          {[0->2], [0->1, 2->1], [1->2, 2->1]}
     *      And so on.
     *  * The class "multiplicative partition holder" keeps all the known partitions,
     *      sorted by the original prime product. The key is a prime decomposition and the
     *      value is a list of multiplicative partitions. For example, for the key [1,1,1]
     *      the values are {{[0->1, 1->1, 2->1]}; {[0->1, 1->1], [2->1]}; {[0->1, 2->1], [1->1]};
     *      {[1->1, 2->1]; [0->1]}; {[0->1], [1->1], 2->1]}}.
     */
    private static class FactorDecomposition    {
        private List<Integer> primePowers;
        public <T> FactorDecomposition(Map<T,Integer> primeDecomposition)    {
            this(primeDecomposition.values());
        }
        public FactorDecomposition(Collection<Integer> primePowers)    {
            this.primePowers=orderDescending(primePowers);
        }
        private List<Integer> orderDescending(Collection<Integer> list)    {
            List<Integer> result=new ArrayList<>(list);
            result.sort(Collections.reverseOrder());
            return result;
        }
        @Override
        public int hashCode()    {
            return primePowers.hashCode();
        }
        @Override
        public boolean equals(Object other)    {
            try    {
                FactorDecomposition pdOther=(FactorDecomposition)other;
                return primePowers.equals(pdOther.primePowers);
            }    catch (ClassCastException cce)    {
                return false;
            }
        }
        public Pair<FactorDecomposition,Integer> oneFactorLess()    {
            List<Integer> parent=new ArrayList<>(primePowers);
            int n=parent.size()-1;
            if (parent.get(n)==1) parent.remove(n);
            else parent.set(n,parent.get(n)-1);
            return new Pair<>(new FactorDecomposition(parent),n);
        }
        @Override
        public String toString()    {
            return primePowers.toString();
        }
        public int totalSize()    {
            int result=0;
            for (int pow:primePowers) result+=pow;
            return result;
        }
        // This shameless exposure shouldn't happen, but eh.
        public List<Integer> getPowers()    {
            return primePowers;
        }
    }

    private static class AbstractProduct implements Comparable<AbstractProduct>    {
        // They key is the index of the prime and the value is the amount of times it appears.
        private SortedMap<Integer,Integer> product;
        public AbstractProduct(SortedMap<Integer,Integer> map)    {
            // Note that this stores the map, doesn't copy it!
            product=map;
        }
        @Override
        public int hashCode()    {
            return product.hashCode();
        }
        @Override
        public boolean equals(Object other)    {
            try    {
                AbstractProduct apOther=(AbstractProduct)other;
                return product.equals(apOther.product);
            }    catch (ClassCastException cce)    {
                return false;
            }
        }
        public static AbstractProduct fromSingleIndex(int index)    {
            return new AbstractProduct(ImmutableSortedMap.of(index,1));
        }
        public AbstractProduct addFactor(int index)    {
            SortedMap<Integer,Integer> newMap=new TreeMap<>(product);
            EulerUtils.increaseCounter(newMap,index);
            return new AbstractProduct(newMap);
        }
        @Override
        public String toString()    {
            return product.toString();
        }
        /*
         * This is necessary so that we have a fixed, predictable ordering. The ordering
         * itself doesn't matter much, but it's later used to guarantee that the "ordered"
         * list has a predictable order, thus preventind duplicate partitions.
         */
        @Override
        public int compareTo(AbstractProduct other) {
            int s=product.size();
            int os=other.product.size();
            if (s!=os) return s-os;
            Iterator<Map.Entry<Integer,Integer>> it1=product.entrySet().iterator();
            Iterator<Map.Entry<Integer,Integer>> it2=other.product.entrySet().iterator();
            for (int i=0;i<s;++i)    {
                Map.Entry<Integer,Integer> e1=it1.next();
                Map.Entry<Integer,Integer> e2=it2.next();
                int dk=e1.getKey()-e2.getKey();
                if (dk!=0) return dk;
                int dv=e1.getValue()-e2.getValue();
                if (dv!=0) return dv;
            }
            return 0;
        }
    }

    private static class MultiplicativePartition    {
        private final List<AbstractProduct> partition;
        public MultiplicativePartition(List<AbstractProduct> products)    {
            // "products" is not copied!
            partition=products;
            partition.sort(null);
        }
        @Override
        public int hashCode()    {
            return partition.hashCode();
        }
        @Override
        public boolean equals(Object other)    {
            try    {
                MultiplicativePartition mpOther=(MultiplicativePartition)other;
                return partition.equals(mpOther.partition);
            }    catch (ClassCastException cce)    {
                return false;
            }
        }
        private MultiplicativePartition addChildAt(int primeIndex,int listIndex)    {
            List<AbstractProduct> listCopy=new ArrayList<>(partition);
            listCopy.set(listIndex,listCopy.get(listIndex).addFactor(primeIndex));
            return new MultiplicativePartition(listCopy);
        }
        private MultiplicativePartition addSingletonPrime(int primeIndex)    {
            List<AbstractProduct> listCopy=new ArrayList<>(partition);
            listCopy.add(AbstractProduct.fromSingleIndex(primeIndex));
            return new MultiplicativePartition(listCopy);
        }
        // This is used to generate new partitions by adding a prime with certain index.
        public void getChildren(int primeIndex,Set<MultiplicativePartition> destination)    {
            for (int i=0;i<partition.size();++i) destination.add(addChildAt(primeIndex,i));
            destination.add(addSingletonPrime(primeIndex));
        }
        @Override
        public String toString()    {
            return partition.toString();
        }
        public FactorDecomposition getAsFactorDecomposition()    {
            Map<AbstractProduct,Integer> counter=new HashMap<>();
            for (AbstractProduct ap:partition) EulerUtils.increaseCounter(counter,ap);
            return new FactorDecomposition(counter);
        }
    }

    private static class MultiplicativePartitionHolder    {
        private SetMultimap<FactorDecomposition,MultiplicativePartition> partitions;
        public MultiplicativePartitionHolder()    {
            partitions=HashMultimap.create();
            // We add the default empty partition for the special case [1] with no primes.
            FactorDecomposition emptyDecomp=new FactorDecomposition(Collections.emptyList());
            MultiplicativePartition emptyPart=new MultiplicativePartition(Collections.emptyList());
            partitions.put(emptyDecomp,emptyPart);
        }
        public Set<MultiplicativePartition> getAllPartitions(FactorDecomposition pd)    {
            if (!partitions.containsKey(pd)) calculatePartitions(pd);
            return partitions.get(pd);
        }
        private void calculatePartitions(FactorDecomposition pd)    {
            Pair<FactorDecomposition,Integer> removedFactor=pd.oneFactorLess();
            Set<MultiplicativePartition> parentPartitions=getAllPartitions(removedFactor.first);
            Set<MultiplicativePartition> newPartitions=new HashSet<>();
            for (MultiplicativePartition oldPartition:parentPartitions) oldPartition.getChildren(removedFactor.second,newPartitions);
            System.out.println("\t"+pd.toString()+" => "+newPartitions.size()+" particiones.");
            partitions.putAll(pd,newPartitions);
        }
        @Override
        public String toString()    {
            return partitions.toString();
        }
        public Map<FactorDecomposition,Integer> getPartitionListAsCountOfAbstractProducts(FactorDecomposition pd)    {
            Map<FactorDecomposition,Integer> result=new HashMap<>();
            for (MultiplicativePartition partition:getAllPartitions(pd))    {
                EulerUtils.increaseCounter(result,partition.getAsFactorDecomposition());
            }
            return result;
        }
        public void clear()    {
            partitions.clear();
        }
    }

    private static class PrimeDecompositionGenerator    {
        private int firstPrimes[];
        public PrimeDecompositionGenerator(int max)    {
            System.out.println("Computing initial prime set...");
            firstPrimes=Primes.firstPrimeSieve(1+max);
        }
        public FactorDecomposition getDecomposition(int in)    {
            Map<Integer,Integer> factors=new HashMap<>();
            while (in>1)    {
                int firstPrime=firstPrimes[in];
                if (firstPrime==0)    {
                    EulerUtils.increaseCounter(factors,in);
                    break;
                }    else    {
                    EulerUtils.increaseCounter(factors,firstPrime);
                    in/=firstPrime;
                }
            }
            return new FactorDecomposition(factors);
        }
        private void clear()    {
            firstPrimes=new int[0];
        }
    }

    private static class FactorialCache    {
        private final Map<Integer,BigInteger> cache;
        public FactorialCache()    {
            cache=new HashMap<>();
            cache.put(0,BigInteger.ONE);
        }
        public BigInteger getFactorial(int k)    {
            BigInteger result=cache.get(k);
            if (result==null)    {
                result=BigInteger.valueOf(k).multiply(getFactorial(k-1));
                cache.put(k,result);
            }
            return result;
        }
    }

    private static class CombinatorialNumberCache    {
        // This stores N*(N-1)*...(N-k) for N fixed and k given.
        private final Map<Integer,BigInteger> cache;
        public CombinatorialNumberCache()    {
            cache=new HashMap<>();
            cache.put(0,BigInteger.ONE);
        }
        public BigInteger getProduct(int k)    {
            BigInteger result=cache.get(k);
            if (result==null)    {
                result=getProduct(k-1).multiply(BigInteger.valueOf(N+1-k));
                cache.put(k,result);
            }
            return result;
        }
    }

    private static Map<FactorDecomposition,Integer> getDecompositionCounter(int n)    {
        PrimeDecompositionGenerator gen=new PrimeDecompositionGenerator(n);
        System.out.println("Generating prime decompositions...");
        Map<FactorDecomposition,Integer> result=new HashMap<>();
        for (int i=1;i<=n;++i) EulerUtils.increaseCounter(result,gen.getDecomposition(i));
        gen.clear();
        return result;
    }

    public static void main(String[] args)    {
        try    {
            PrintStream ps=new PrintStream("C:\\out452_3.txt");
            System.setOut(ps);
        }    catch (Exception exc)    {
            System.out.println("D'oh!");
        }
        // First we create an object which counts how many times does each distinct prime
        // decomposition appear.
        Map<FactorDecomposition,Integer> decompCounter=getDecompositionCounter(N);
        // Then, using the capabilities of the MultiplicativePartitionHolder, we translate
        // decompCounter into an object which counts how many partitions have certain amount
        // of elements.
        System.out.println("Sorting partitions by size...");
        MultiplicativePartitionHolder partitionHolder=new MultiplicativePartitionHolder();
        Map<FactorDecomposition,Long> partitionCounter=new HashMap<>();
        for (Map.Entry<FactorDecomposition,Integer> entry:decompCounter.entrySet())    {
            long multiplier=entry.getValue();
            Map<FactorDecomposition,Integer> partitionsBySize=partitionHolder.getPartitionListAsCountOfAbstractProducts(entry.getKey());
            for (Map.Entry<FactorDecomposition,Integer> howManyPartitions:partitionsBySize.entrySet()) EulerUtils.increaseCounter(partitionCounter,howManyPartitions.getKey(),multiplier*howManyPartitions.getValue());
        }
        decompCounter.clear();
        partitionHolder.clear();
        System.out.println("Performing the last sum...");
        FactorialCache factorials=new FactorialCache();
        CombinatorialNumberCache combinatorials=new CombinatorialNumberCache();
        BigInteger sum=BigInteger.ZERO;
        for (Map.Entry<FactorDecomposition,Long> entry:partitionCounter.entrySet())    {
            BigInteger multiplier=BigInteger.valueOf(entry.getValue());
            FactorDecomposition decomp=entry.getKey();
            BigInteger augend=combinatorials.getProduct(decomp.totalSize());
            for (int pow:decomp.getPowers()) augend=augend.divide(factorials.getFactorial(pow));
            augend=augend.multiply(multiplier);
            sum=sum.add(augend);
        }
        sum=sum.mod(BigInteger.valueOf(MOD));
        System.out.println(sum);
    }
}
//projecteuler=452
