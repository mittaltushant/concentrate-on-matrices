---
title: Matrix Chernoff Bounds for sums of independent random variables 
date: 2023-01-26
category: notes
speaker: Tushant Mittal
abstract: "We do a matrix-version of our dear Chernoff bound."
---


## Summary
In this talk, we proved a matrix-version of a simple Chernoff bound. The scalar version which we will generalize is the following. 
### Scalar Chernoff
Let $X = \sum_i \eps_ia_i $ where $a_i \in \R$ are fixed, and, $\eps_i$ are $\pm 1$ Rademacher random variables. Then,
\\[ \prob{\sum_i \eps_ia_i > t} ~\leq~ e^{\frac{-t^2}{2\sigma^2}},\\;\\; \sigma^2 = \ex{X^2} =\\; \sum_i a_i^2
\\]
Let's look at a standard proof of this to see where it breaks down. For a single Rademacher variable, one can show that $\ex{e^{\theta \eps a}} = \mathrm{cosh}(\theta a) \leq e^{\frac{\theta^2 a^2}{2}}$. The goal is to break the sum on $n$ such variables as a product of these. Fix $\theta > 0\in \R$.
\\[
\begin{align}
\prob{X > t} ~&=~ \prob{e^{\theta X} > e^{\theta t}}\\;\\;\small{\text{ (Monotonicity of exponential)} }\\\
 &\leq e^{-\theta t}\\,\ex{ e^{\theta X}} \\;\\;\small{\text{ (Markov)} }\\\
 ~&=~ e^{-\theta t}\\,\ex{ \prod_i e^{\theta \eps_i a_i}} \\;\\; \small{(e^{a+b} = e^a\cdot e^b)} \\\
 ~&=~ e^{-\theta t}\\,\prod_i\ex{ e^{\theta \eps_i a_i}} \\;\\; \small{\text{ (Independence)} } \\\
 ~&\leq~ e^{-\theta t}\\,\prod_i e^{\frac{\theta^2 a_i^2}{2}}  \\;\\; \small{\text{ (Plugging in the 1-variable bound )} }\\\
 ~&\leq~ e^{\frac{-t^2}{2(\sum_i a_i^2)}}  \\;\\; \small{\text{ (Optimizing for }\theta) }
\end{align}
\\]

### Matrix Chernoff
We will generalize the above setup as follows.
Let $X = \sum_i \eps_iA_i $ where $A_i$ are fixed Hermitian matrices, and, $\eps_i$ are $\pm 1$ Rademacher random variables. Then,
\\[ \prob{\lambda_{\max}(X) > t} ~\leq~ e^{\frac{-t^2}{2\sigma^2}} 
\\]

We will give two proofs each of which give a slightly different $\sigma$. The first question is to make sense of exponentials of matrix valued random variables. This, fortunately, is easy. 

### Lifting functions to matrices
Let $f: I \to \R$ be a function defined on an interval $I$. Let $A = Q\Lambda Q^* $ be an Hermitian matrix such that all eigenvalues of $A$ lies in $I$, i.e., $\mathrm{Spec}(A)\subseteq I$. Then, we can define $f(A) := Qf(\Lambda) Q^*$, where $f(\Lambda)$ is obtained by applying $f$ entry-wise to the diagonal matrix $\Lambda$, i.e., $f(\Lambda) = \mathrm{diag}(f(\lambda_1), \cdots, f(\lambda_d))$.

Thus, we can talk of $e^A$ for any matrix $A$, of $\log (B)$ for any positive-definite matrix $B$. In particular, $\log (\sum_i \alpha_i e^{B_i})$ is well-defined for $\alpha_i > 0$ as $e^{B_i}$ are all positive-definite and the set of PD matrices forms a (open) cone.   
### Mimicking the scalar proof

 One can repeat the scalar argument for the first two steps but it is not clear how to handle the term, $\ex{\lambda_{\max}(e^{\theta X})}$. Ideally, we would like to have $\prod_i \ex{\lambda_\max (e^{\theta \eps_i A_i})}$. The key difficulty is that $e^{A+B}\neq e^Ae^B$ for matrices. This can be resolved (by at least) two approaches. 

\\[
\begin{align}
\small{[\text{Naive Hope}]}\\;\\;\\; &\ex{\lambda_{\max}(e^{\theta X})} ~=~ \ex{\lambda_{\max}\left(\prod_i e^{\theta \eps_i A_i}\right)} ~\leq~ \prod_i \ex{\lambda_\max (e^{\theta \eps_i A_i})}  \\\
\small{[\text{Ahlswede-Winter}]}\\;\\;\\; &\ex{\lambda_{\max}(e^{\theta X})} ~\leq~ \ex{\tr(e^{\theta X})} ~\leq~  d\prod_i \ex{\lambda_\max (e^{\theta \eps_i A_i})}\\\
\small{[\text{Tropp}]}\\;\\;\\; &\ex{\lambda_{\max}(e^{\theta X})} ~\leq~ \ex{\tr(e^{\theta X})} ~\leq~   \tr\\, \mathrm{exp}\left( \sum_i \log \ex{(e^{\theta \eps_i A_i})}\right) \\\
\end{align}
\\]

### Using the trace inequalities

- [Ahlswede-Winter] By definition of the matrix exponential, $\lambda_\max(e^A) = e^{\lambda_\max(A)}$. Thus, treating $\lambda_\max(A)$ as a scalar random variable, we can plug it into the scalar inequality we used earlier, $\ex{e^{\theta \eps a}} \leq e^{\frac{\theta^2 a^2}{2}}$. Thus, we get, $\ex{\lambda_\max (e^{\theta \eps_i A_i})} \leq e^{\frac{\theta^2 \lambda_\max(A_i)^2}{2}}$. It is now exactly like the scalar Chernoff and we get a variance term of $\sum_i \lambda_\max(A_i)^2$.

- [Tropp] We need two more facts. 
    - Firstly, a matrix version of the 1-variable inequality. This is given by $\log \ex{e^\theta \eps_i A_i} \preceq \frac{\theta^2 A_i^2}{2}$. Here the order being used is the Loewner order ($A\preceq B$ iff $B-A$ is PSD). The proof is analogous to the scalar proof and is given in Tropp's book [].
    -  The fact that trace-exponential is monotone in the sense that if $A\preceq B$, $\tr\\, e^A \leq \tr\\, e^B$. This is easy to establish by using the fact that $A\preceq B$ implies $\lambda_i(A)\leq \lambda_i(B)$ which itself follows by Courant-Fischer. 
    
Now, we are ready. 
\\[
\begin{align}
\log \ex{e^\theta \eps_i A_i} &~\preceq~ \frac{\theta^2 A_i^2}{2}\\;\\;\\; \small{\text{[Fact 1]}}\\\
\sum_i \log \ex{e^\theta \eps_i A_i} &~\preceq~ \frac{\theta^2 \sum_i A_i^2}{2}\\\
\tr \\; \mathrm{exp}\left(\sum_i \log \ex{e^\theta \eps_i A_i}\right) &~\leq~ \tr \\; \mathrm{exp}\left(\frac{\theta^2 \sum_i A_i^2}{2}\right) \\;\\; \\; \small{\text{[Fact 2]}}\\\
&~\leq~ d\\, \lambda_\max \left( \mathrm{exp}\left( \frac{\theta^2 \sum_i A_i^2}{2}\right)\right)\\\
&~=~ d\\, \mathrm{exp}\left( \frac{\theta^2 \lambda_\max\left(\sum_i A_i^2\right)}{2}\right).
\end{align}
\\]
This, gives us a variance term of $\lambda_\max\left(\sum_i A_i^2\right) $ which is better than the earlier one of $\sum_i \lambda_\max\left(A_i\right)^2$ by up to a factor of $d$. This matters as we have the factor of $d$ in the exponent.

### Deriving the trace inequalities
- [Golden-Thompson] Ahlswede--Winter use the following sequence of inequalities, 
	\\[\tr(e^{A+B}) \leq \tr(e^A e^B) \leq \lambda_\max(e^A)\tr(e^B).\\]
	
	The first inequality is called the **Golden-Thompson** inequality which can be derived from the Lie-Trotter formula which says that although $e^{A+B}\neq e^Ae^B$, this holds in the limit. Formally,  $e^{A+B} = \lim_{n\to \infty} \\left(e^{A/n}e^{B/n} \\right)^{n}$. However, the Golden-Thompson inequality is false for three or more matrices. Thus, it cannot be applied directly. To overcome this, Ahlswede--Winter use the second inequality, a proof of which can be found in Harvey's [notes](https://www.math.uwaterloo.ca/~harvey/W11/Lecture11Notes.pdf). We now show the derivation, 
\\[ 
\begin{align}
\ex{\tr(e^{\theta X + 0)})} ~=~ \ex{\tr(e^{\theta \sum_i \eps_i A_i + 0})} ~&\leq~ \ex{\lambda_\max(e^{\eps_1A_1}) \tr(e^{\theta \sum_{i=2}^n \eps_i A_i+ 0})}\\\
 ~&\leq~ \ex{\prod_{i} \lambda_\max(e^{\eps_iA_i}) \tr(e^{0})} \\\
~&\leq~ d \prod_i \ex{ \lambda_\max(e^{\eps_iA_i})}
\end{align}
\\]

- [Lieb's Concavity] Tropp's insight is that one must instead work with the cumulant generating function, $\log e^{X}$. The advantage of this POV is that this approach generalizes to a much more general settings. Moreover, as we have seen, it gives tighter bounds. The key is the following result of Lieb. 

{{< css.inline >}}
<span class="theorem">(Lieb) <i> Let $H$ be any Hermitian matrix. Then, the function $f(A) = \tr\\, \exp (H +\log A)$ defined on the cone of
positive-definite matrices, is concave. </i>
</span>
{{</ css.inline >}}

One now only needs to apply Jensen's to get the trace inequality, 
\\[
\begin{align}
\mathbb{E}\Ex{X_n}{\tr\\, e^{\sum_{i<n} X_i+ X_n} } ~&=~ \mathbb{E}\Ex{X_n}{\tr\\, e^{\sum_{i<n} X_i+ \log (e^{X_n})} }\\\
~&\leq~ \mathbb{E}\\, \tr\\, e^{\sum_{i<n} X_i +\log \Ex{X_n}{e^{X_n}}} \\;\\;\small{\text{(Jensen's)}}\\\
~&\leq~ \tr\\, e^{\sum_{i} \log \Ex{X_i}{e^{X_i}}} \\;\\;\small{\text{(Repeating the above steps for all variables)}}
\end{align}
\\]  
 

## Resources 
- Tropp, Joel A. "An Introduction to Matrix Concentration Inequalities." [arXiv](https://doi.org/10.48550/arXiv.1501.01571). 
- [Lecture Notes](https://www.math.uwaterloo.ca/~harvey/W11/Lecture11Notes.pdf) on Ahlswede-Winter Inequality by Nicholas Harvey.
  
- Garg, Ankit, Lee, Yin T., Song, Zhao, and Nikhil Srivastava. "A Matrix Expander Chernoff Bound." [arXiv](https://doi.org/10.48550/arXiv.1704.03864).

- Talk by Joel Tropp at [[Youtube Link]](https://www.youtube.com/watch?v=T9ViSznHeUE)