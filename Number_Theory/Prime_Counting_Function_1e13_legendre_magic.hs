-- This file is a Legendre prime counting function...
-- credits : GordonBGood
-- tranrslated from cpp source: misterkeemu
-- Tested on: Yosupo Counting Primes
-- Complexity: O(n^(3/4)/((log n)^2))
-- Returns number of primes <= n

{-# OPTIONS_GHC -O2 #-}
{-# LANGUAGE FlexibleContexts, BangPatterns #-}

import Data.Int ( Int64, Int32 )
import Data.Bits ( shiftL, shiftR, (.&.), (.|.) )
import Control.Monad ( forM_, when )
import Control.Monad.ST (ST, runST)
import Data.Array.Base ( STUArray(..), unsafeAt,
                         castSTUArray, unsafeFreezeSTUArray,
                         MArray(unsafeNewArray_, unsafeRead, unsafeWrite) )

primeCount :: Int64 -> Int64
primeCount n =
  if n < 3 then (if n < 2 then 0 else 1) else
  let
    {-# INLINE divide #-}
    divide :: Int64 -> Int64 -> Int
    divide nm d = truncate $ (fromIntegral nm :: Double) / fromIntegral d
    {-# INLINE half #-}
    half :: Int -> Int
    half x = (x - 1) `shiftR` 1
    rtlmt = floor $ sqrt (fromIntegral n :: Double) 
    mxndx = (rtlmt - 1) `div` 2
    (!nbps, !nrs, !smalls, !roughs, !larges) = runST $ do
      -- becomes `smalls` LUT -> the current counts of odd primes to index...
      mss <- unsafeNewArray_ (0, mxndx) :: ST s (STUArray s Int Int32)
      let msscst =
            castSTUArray :: STUArray s Int Int32 -> ST s (STUArray s Int Int64)
      mdss <- msscst mss -- for use in adjing counts LUT
      forM_ [ 0 .. mxndx ] $ \ i -> unsafeWrite mss i (fromIntegral i)
      -- becomes `roughs` LUT -> the current "k-roughs" for base prime sieved...
      mrs <- unsafeNewArray_ (0, mxndx) :: ST s (STUArray s Int Int32)
      forM_ [ 0 .. mxndx ] $ \ i -> unsafeWrite mrs i (fromIntegral i * 2 + 1)
      -- becomes `larges` LUT -> the current count of odd primes indexed for
      -- the inverse of the current "k-roughs" in the table above...
      mls <- unsafeNewArray_ (0, mxndx) :: ST s (STUArray s Int Int64)
      forM_ [ 0 .. mxndx ] $ \ i ->
        let d = fromIntegral (i + i + 1)
        in unsafeWrite mls i (fromIntegral (divide n d - 1) `div` 2)
      cmpsts <- unsafeNewArray_ (0, mxndx) :: ST s (STUArray s Int Bool)
      -- partial sieves to quad root of counting range, adjusting and
      -- accumulating LUT's so that the overall current results are
      -- accumulated to the `mls`/`larges` array...
      -- also outputs `cbpi`/`nbps` is the number of base prime sieved and
      -- `rlmti`/`nrs` is the effective size of the "k-roughs" sized LUT's...
      let loop i !cbpi !rlmti =
            let sqri = (i + i) * (i + 1) in
            if sqri > mxndx then do
              fss <- unsafeFreezeSTUArray mss
              frs <- unsafeFreezeSTUArray mrs
              fls <- unsafeFreezeSTUArray mls
              return (cbpi, rlmti + 1, fss, frs, fls)
            else do
              v <- unsafeRead cmpsts i
              if v then loop (i + 1) cbpi rlmti else do
                unsafeWrite cmpsts i True -- cull current bp so not a "k-rough"!
                let bp = i + i + 1
                    -- partial cull by current base prime...
                    cull c = if c > mxndx then return () else do
                               unsafeWrite cmpsts c True; cull (c + bp)
                    -- adjust `mls` array for current partial sieve;
                    -- also adjusts effective sizes of `mrs` and `mls`...
                    part ri nri = -- old "rough" index to new one...
                      if ri > rlmti then return (nri - 1) else do
                        r <- unsafeRead mrs ri -- "rough" always odd!
                        t <- unsafeRead cmpsts (fromIntegral r `shiftR` 1)
                        if t then part (ri + 1) nri else do -- skip newly culled
                          olv <- unsafeRead mls ri
                          let m = fromIntegral r * fromIntegral bp
                          -- split -> when multiple <= square root:
                          --   quotient `n / m` will be less than square root so
                          --   `mls` index will be found from indexing `mss`
                          --   (adjusted by current number bp's not in `mls`)...
                          adjv <- if m <= fromIntegral rtlmt then do
                                    let ndx = fromIntegral m `shiftR` 1
                                    sv <- unsafeRead mss ndx
                                    unsafeRead mls (fromIntegral sv - cbpi)
                          -- else quotient will be less than square root so
                          --   quotient can be directly indexed from `mss`...
                                  else do
                                    sv <- unsafeRead mss (half (divide n m))
                                    return (fromIntegral sv)
                          -- move "rough" and new "large" values to new places:
                          -- adjv includes number base primes already in `olv`
                          unsafeWrite mls nri (olv - (adjv - fromIntegral cbpi))
                          unsafeWrite mrs nri r; part (ri + 1) (nri + 1)
                    !pm0 = ((rtlmt `div` bp) - 1) .|. 1 -- max base prime mult
                    -- adjust `mss` counting table for current partial sieve;
                    -- for array range to `lmti`; prime multiple to `pm`...
                    -- adjust 64-bits at a time where possible for speed...
                    adjc lmti pm =
                      if pm < bp then return () else do
                        c <- unsafeRead mss (pm `shiftR` 1)
                        let ac = c - fromIntegral cbpi -- correction
                            bi = (pm * bp) `shiftR` 1 -- start array index
                            adj si = if si > lmti then adjc (bi - 1) (pm - 2)
                                     else do ov <- unsafeRead mss si
                                             unsafeWrite mss si (ov - ac)
                                             adj (si + 1)
                            ac64 = fromIntegral ac :: Int64
                            dac = (ac64 `shiftL` 32) .|. ac64
                            dbi = (bi + 1) `shiftR` 1
                            dlmti = (lmti - 1) `shiftR` 1
                            dadj dsi = if dsi > dlmti then return ()
                                      else do dov <- unsafeRead mdss dsi
                                              unsafeWrite mdss dsi (dov - dac)
                                              dadj (dsi + 1)
                        when (bi .&. 1 /= 0) $ do
                          ov <- unsafeRead mss bi
                          unsafeWrite mss bi (ov - ac)
                        dadj dbi
                        when (lmti .&. 1 == 0) $ do
                          ov <- unsafeRead mss lmti
                          unsafeWrite mss lmti (ov - ac)
                        adjc (bi - 1) (pm - 2)
                cull sqri; nrlmti <- part 0 0; adjc mxndx pm0
                loop (i + 1) (cbpi + 1) nrlmti
      loop 1 0 mxndx
    !ans0 = unsafeAt larges 0 - -- combine all counts; each includes nbps...
              sum [ unsafeAt larges i | i <- [ 1 .. nrs - 1 ] ]
    -- adjust for all the base prime counts subracted above...
    !adj = (nrs + 2 * (nbps - 1)) * (nrs - 1) `div` 2
    !adjans0 = ans0 + fromIntegral adj
    -- add counts for base primes above quad root counting range
    -- to cube root counting range multiplied by rough primes above
    -- the base prime as long as the quotient of `n` divided by the
    -- multiple is greater than the base prime; counts of indexed by
    -- the quotient as above...
    -- since all `roughs` are now prime, the multiple will always be
    -- just two primes so the compensation will always be added;
    -- also, the product will always be > the square root of the range so
    -- the quotient will always be less than the square root of the range and
    -- only the `smalls` count LUT needs be used (second case from above loop).
    loopr ri !acc =
      if ri >= nrs then acc else
      let r = fromIntegral (unsafeAt roughs ri)
          q = n `div` r
          lmtsi = half (fromIntegral (q `div` r))
          lmti = fromIntegral (unsafeAt smalls lmtsi) - nbps
          addcnt pi !ac =
            if pi > lmti then ac else
            let p = fromIntegral (unsafeAt roughs pi)
                ci = half (fromIntegral (divide q p))
            in addcnt (pi + 1) (ac + fromIntegral (unsafeAt smalls ci))
      in if lmti <= ri then acc else
         -- adjust for the `nbps`'s over added in the `smalls` counts...
         let !adj = fromIntegral ((lmti - ri) * (nbps + ri - 1))
         in loopr (ri + 1) (addcnt (ri + 1) acc - adj)
  in loopr 1 adjans0 + 1 -- add one for only even prime of two!

main :: IO ()
main = print . primeCount . read =<< getLine
