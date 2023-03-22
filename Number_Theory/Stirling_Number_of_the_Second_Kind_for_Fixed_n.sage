for a fixed n and all 0≤𝑘≤n is to use the generating function identity ∑n𝑘=0𝑆(n,𝑘)𝑡𝑘⎯⎯=𝑡n. Here 𝑡𝑘⎯⎯ is a falling factorial.
SageMath's built-in stirling_number2 method uses Maxima or GAP's implementations;
Maxima seems faster in a little testing. You're not going to get to n in the "millions". Already for n=2000 the largest coefficient has 4348 digits, well past 32-bit.