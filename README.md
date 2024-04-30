### A modification of Balas's Additive Algorithm, specialized for the following binary optimization problem:

<br>

&nbsp;&nbsp;&nbsp;&nbsp;The goal is to pick the minimum number of objects - where each object has
an assigned score in three categories - such that the sum of the scores across
all three categories for the selected objects is greater than some threshold, $δ$.
The selected group of objects is denoted as the set $K$, where each object in
the group has some defined score in three categories: $k_a,k_b,k_c ∈ [0,1]$.

$$
min\ |K|
$$

$$subject\ to:$$

$$
\sum_{k∈K} k_i ≥ δ, i ∈ {a,b,c}
$$
