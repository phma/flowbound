# A336433
This program computes the sequence of the number of sequences of *n* numbers from 1 to *n* that do not have a subsequence that adds up to *n*.

I thought of this sequence when planning how to load a GPU with blocks of data points to crunch. Each task block has a multiple of 1024 points (the ones that don't are processed by the CPU) and needs to be crunched down to a few numbers. My GPU has 36 groups of 256 cores each, and I could write the kernel so that each core processes 4 points, and then the cores in each group add up the results.

I'd prefer to break task blocks as little as possible, while keeping the GPU as full as possible. (The name of the program is short for "full processor".) So if a task block is 87 kibipoints, I'll break it into two blocks of 36 kibipoints and one of 15 kibipoints, which I can process together with a 21-kibipoint piece of another task block.

Given some random task block sizes, what is the probability that one cannot fill the processor with *n* work groups, assuming that remainders (15 above) are processed first? That is what this sequence answers, when divided by *n*^*n*. This problem is theoretical, as finding a subsequence of *n* numbers from 1 to *n* that adds up to *n* is a hard problem, and the computer has to feed the GPU with data possibly thousands of times a second. But finding two of the *n* numbers that add up to *n* is easy, and most of the time, there is such a pair.