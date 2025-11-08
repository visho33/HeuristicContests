# Solution to [HOT 2024 Round 1](https://algotester.com/en/Contest/Display/110664)

## Contest Placement
I participated in the contest and ranked **6th place**. This repository contains my solution, which I hope can serve as a useful reference.

## Approach

The core idea of my solution is to implement a **score function** that tracks the exact state of the pipeline system. This allows me to evaluate the score precisely, rather than relying on a greedy approximation. Here's an overview of my approach:

### Key Concepts and Steps

1. **Score Function:**
   I implemented the score function to efficiently update the score and ensure accurate tracking when changes occur in the pipeline. The function supports adding and removing requests, as well as virtual operations for faster processing.

2. **Initial Insights:**
   - The first insight was that it's always optimal to set `b = 1` in the score function, which simplifies the calculation.
   - The second insight was that maximizing the number of pipelines is generally advantageous, so I analyzed the formulas to support this.

3. **Greedy Approach:**
   Requests are added one by one using a greedy approach, where each request is assigned to the pipeline that minimizes the score function. Initially, I set `tau = 0` to avoid having `tau` dominate over the other components of the latency function, ensuring a smooth start.

4. **Optimizing Pipeline Count:**
   After setting up the initial solution, I worked to reduce the number of pipelines. I did this by checking if fewer pipelines resulted in a better score, helping to optimize the pipeline system.

5. **Simulated Annealing:**
   To further improve the solution, I employed **simulated annealing** to make local changes between pipelines. I swapped requests between pairs `(j1, r1)` and `(j2, r2)`, first trying 45 swaps with fixed jobs and random requests, and then 45 random swaps.
   - Every 90 local changes, I performed a **global change**, where I erased one request from each of the pipelines with the highest latency and reassigned those requests greedily.

6. **Efficiency Considerations:**
   I ensured that the score was updated only when necessary. This efficient approach minimized redundant calculations, particularly when small changes (like swapping requests) did not significantly affect the overall score.

7. **Improvements and Future Work:**
   Although the solution is solid, I believe that optimizing the parameters for the simulated annealing process could improve the results further. Unfortunately, due to time constraints during the competition, I couldn't fine-tune these parameters. However, I believe that this could help improve performance and lead to even better results.

## Solution Code

The code provided in this repository implements the approach described above.
