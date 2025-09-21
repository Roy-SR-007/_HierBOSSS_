# _HierBOSSS_: Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions

<p align="center">
  <img src="hierbosss_tree.gif" alt="HierBOSSS_logo" width="800"/>
</p>

This repository holds the source code and implementation of **HierBOSSS** for Bayesian structural learning of scientific symbolic expressions.

<br><br>

<figure align="center">
  <img src="assets/symbolic_tree_representation.png" alt="symbolic_tree_representation" width="800"/>
  <figcaption><em>Figure 1: Symbolic tree representation of scientific expressions.</em></figcaption>
</figure>

<br><br>


**HierBOSSS** models symbolic expressions through an operator-induced sum-of-symbolic trees. Conjugate priors are assigned to model regression parameters, while a regularizing prior is designed for the individual symbolic tree structures. To perform inference from the HierBOSSS-induced posterior distribution, we develop an efficient Metropolis-within-partially-collapsed Gibbs Markov chain Monte Carlo (MCMC) sampling algorithm. This GitHub repository showcases the success of HierBOSSS in discovering interpretable scientific laws. Specifically, we demonstrate HierBOSSS' ability to recover and learn well-known physics-based Feynman equations and identify meaningful descriptors in single-atom catalysis. Moreover, HierBOSSS consistently outperforms state-of-the-art symbolic regression methods (which often suffer from training data size, noise, and overly complicated symbolic expressions), offering advantages both in symbolic expression discovery and computational efficiency.

<br><br>

<figure align="center">
  <img src="assets/HierBOSSS_comparison.png" alt="HierBOSSS_comparison" width="800"/>
  <figcaption><em>Figure 2: Comparing HierBOSSS with competing symbolic regression methods in learning the Coulomb's law.</em></figcaption>
</figure>
