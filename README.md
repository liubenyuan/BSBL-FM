# BSBL-FM #

This is a fast implementation of the Block Sparse Bayesian Learning [BSBL](https://sites.google.com/site/researchbyzhang/) algorithm. The developed algorith is based on the Fast Marginalized **(FM)** likelihood maximization algorithm, which yields ~8 times speedup while also pertains nearly the same recovery performances.

# A Short Introduction #

A CS algorithm aims to solve, **Y** = **Phi** **X** + **N**, where **Y** is the measurement matrix of size M times T, **Phi** is the under-determined sensing matrix of size M times N, **X** is the signal.

Compressed sensing, can recover **X** given **Y** and the under-determined matrix **Phi**. When **T=1**, we called it Single Measurement Vector (**SMV**) Model, with **T>1**, it is the Multiple Measurement Vector (**MMV**) model.

Block Sparse, assumes that **x** can be partitioned into blocks **x** = {**x_1**, ... , **x_g**}. The non-zero entries cluster within some blocks and *zeros* otherwise. If *d* out of *g* blocks are non-zero, then the block sparsity is defined as *d/g*. Exploiting both the block sparsity and the intra-block correlation is the source of the magic of all BSBL algorithms.

Our **BSBL-FM** algorithm, is an ultra fast implementation of the original BSBL framework, which brings about **~8** times speedup. What's more, It can worked in all the scenarios include:

> -    SMV sparse
> -    MMV sparse
> -    SMV block sparse
> -    MMV block sparse
> -    Real-valued
> -    Complex-valued

See the demos and implementations below for more details.

# Codes and Data #

The `.m` codes are:

> **CODE:**
> 
> - **BSBL_FM.m**: the main algorithm, also called **MBSBL-FM** in MMV model
> - **BSBL_BO.m**: Zhilin's BSBL-BO algorithm.
> - **demo_smv_real.m**: the real valued SMV block sparse demo
> - **demo_smv_complex.m**: the complexed valued SMV block sparse demo
> - **demo_mmv.m**: the real valued MMV block sparse demo
> - **demo_fecg.m**: the demo code for FECG dataset recovery

The `.mat` data files are:

> **DATA:**
>
> - **demo.mat**: the data for SMV case, contains *re*, *im* vectors
> - **signal_01.mat**:  FECG datasets used in BSBL-BO by Zhilin
> - **Phi.mat**:  the sensing matrix for CS FECG data

# Citations #

If you find the **BSBL-FM** algorithm useful, please cite:

```bibtex
@Article{liu2013energy,
    Title = {Energy Efficient Telemonitoring of Physiological Signals via Compressed Sensing: A Fast Algorithm and Power 
Consumption Evaluation},
    Author = {Liu, Benyuan and Zhang, Zhilin and Xu, Gary and Fan, Hongqi and Fu, Qiang},
    Journal = {Biomedical Signal Processing and Control},
    Year = {2014},
    Pages = {80--88},
    Volume = {11C}
}
```

```bibtex
@InProceedings{liu2013compression,
    Title = {Compression via Compressive Sensing: A Low-Power Framework for the Telemonitoring of Multi-Channel Physiological 
Signals},
    Author = {Benyuan Liu and Zhilin Zhang and Hongqi Fan and Qiang Fu},
    Booktitle = {2013 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)},
    Year = {2013},
    Organization = {IEEE},
    Pages = {9--12}
}
```

More powerful **STSBL** algorithm developed by Zhilin Zhang is available on-line:
```bibtex
@InProceedings{ZhangAsilomar2013,
    Title = {Compressed Sensing for Energy-Efficient Wireless Telemonitoring: Challenges and Opportunities},
    Author = {Zhilin Zhang and Bhaskar D. Rao and Tzyy-Ping Jung},
    Booktitle = {Asilomar Conference on Signals, Systems, and Computers (Asilomar 2013)},
    Year = {2013}
}
```

```bibtex
@Article{zhang2014spatiotemporal,
    Title = {Spatiotemporal Sparse Bayesian Learning with Applications to Compressed Sensing of Multichannel EEG for Wireless 
Telemonitoring and Brain-Computer Interfaces},
    Author = {Zhilin Zhang and Tzyy-Ping Jung and Scott Makeig and Bhaskar D. Rao and Zhouyue Pi},
    Journal = {(Accepted) IEEE Trans. on Neural Systems and Rehabilitation Engineering},
    Year = {2014}
}
```

