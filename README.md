# MITM Attack on the 7-round Whirlpool (small version with state size 128)

This repository provides the implementation of the MITM attack on the 7-round Whirlpool 
(small version with 4 x 4 state, of size 128 bits).

Round functions of AES are used to simulate both the round functions and key-schedule of 
small version of Whirlpool compression function (round constants are from full Whirlpool).

## On verification of the whole attack:

* To verify the computational complexity of the attack, we count the number of calls of the forward computations, backward computations, and full computations of 7-round hashing. They are denoted with prefix `complexity_F`, `complexity_B`, and `complexity_M` in the program. The total complexity is computed as `(complexity_F + complexity_B)/2 + complexity_M`.

* To make the verification feasible in limited time, we find matches on `PARTIAL` bits instead of 128 bits. In the program, one can set `PARTIAL` = `36`, `40`, `44`, `48`... The concrete patterns of the partial matches are defined using `PARTIAL_MATCH_TARGET_MASK` in [defines.h] (we give priority to the **first diagonal** and the **first column** of the state).


## On precomputation:

* We emphasize that even if the precomputation is included in the main attack procedure, the computational 
complexity will stay at the same order of magnitude. 

* The degree of freedom for the forward is five bytes (2^40). However, we use four bytes (2^32) only, because the degree of freedom for backward is the bottleneck. 

* In the program, `DoFF` is the degree of freedom for the forward in bits. Setting `DoFF` to be 32, 
the computational complexity of the precomputation for the forward chunk is 2^32. 
In our implementation, the memory requirement of precomputation is about 64 GB. 
Together with memory used in the attack, each thread can use up to 256 GB.

* To test the program on a PC with less than 256 GB memory, one can set `DoFF` to be 28. Then, the memory requirement is reduced (but still not small), while the computational complexity for the whole attack will be increased by a factor of 2^4. For example,  when setting `DoFF` to be 32, the expected computational complexity of finding matches on 40-bit is 2^32, while the complexity will be 2^36 if one set `DoFF` to be 28.


## On the results:

* The attacks were verified on a server with 64 cores (use only 1 core for 36~44 bits matching and 4 cores for 48 bits matching) and 1TB memory. 
  * To find a match on 36 bits, it took 29.3653 mins (see [F36.log]);
  * To find a match on 40 bits, it took 32.7435 mins (see [F40.log]);
  * To find a match on 44 bits, it took 518.294 mins (see [F44.log]);
  * To find a match on 48 bits, it took 4979.08 mins (see [F48.log]).
  
  It used 4 OpenMP threads in parallel for matching on 48 bits. 
  The computational complexity is the **sum** of that counted by each thread, 
  and the time is the **combined CPU time of all threads** as returned by clock().


## References

[1] Zhenzhen Bao, Jian Guo, Danping Shi, and Yi Tu: Superposition Meet-in-the-Middle Attacks: Updates on Fundamental Security of AES-like Hashing. Cryptology ePrint Archive, Report 2021/575, 2021. https://eprint.iacr.org/2021/575
