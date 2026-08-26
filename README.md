# M-high-learner

The code for analyzing simulation and experimental data is located in the `inst` folder.

Please install R packages "ACAT" firstly.

ACAT package link:https://github.com/yaowuliu/ACAT


## 🎥 Video Presentation

<div align="center">
  <a href="https://youtu.be/0cCmMTm6Cj8">
    <img src="https://img.youtube.com/vi/0cCmMTm6Cj8/maxresdefault.jpg" alt="M-High-Learner Video Presentation" width="85%">
  </a>
</div>

<p align="center">
  🍿 <b>Click the image above to watch the full video presentation on YouTube.</b>
</p>

<p align="center">
  ⏱️ <b>Video Timestamps:</b>
  <br>
  <code>00:00</code> Motivation & Background &nbsp;|&nbsp;
  <code>00:44</code> Problem Formulation & Causal Framework &nbsp;|&nbsp;
  <code>01:48</code> Three-Stage M-High-Learner Pipeline &nbsp;|&nbsp;
  <code>03:22</code> Subgroup Clustering & Decision Rules &nbsp;|&nbsp;
  <code>04:20</code> Real-Data Application: Framingham Heart Study &nbsp;|&nbsp;
  <code>04:54</code> Conclusion & Future Work
</p>

---

### 💡 Overview & Key Highlights

* **Framework:** **M-High-Learner** is designed to discover heterogeneous mediation effects when mediators are high-dimensional omics measurements by estimating Conditional Average Indirect Effects (CAIE).
* **Core Pipeline:**
  1. **Filtering & Screening:** $R^2_{\text{med}}$-based screening followed by a single-mediator M-Learner to select relevant candidate mediators.
  2. **Joint Prediction:** Multivariate machine learning models for vector-valued mediators to account for co-regulatory networks.
  3. **Subgroup Discovery:** Distance-based clustering (modified t-SNE + C-means) paired with interpretable decision trees to derive clinical subgroup rules.
* **Empirical Validation:** Demonstrated on the *Framingham Heart Study*, identifying 4 distinct BMI- and age-defined subgroups with strong mediation signals ($R^2_{\text{med}} = 0.057, \text{SOS} = 0.36$).

---

### 🔗 Resources & Citation

- 📄 **Paper:** [Link to Paper/ArXiv](https://arxiv.org)
- 🎬 **Video Presentation:** [YouTube Link](https://youtu.be/0cCmMTm6Cj8)
- 📊 **Presentation Slides:** [Download PDF](assets/slides.pdf)










If you are interested in the method and package or report the bugs, please cantact with Xingyu Li: lxingyu1996@gmail.com
