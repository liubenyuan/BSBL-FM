The demo codes for the
    << Fast Marginalized Block Sparse Bayesian Method >>
can work with

    SMV sparse
    MMV sparse
    SMV block sparse
    MMV block sparse
    Real-valued
    Complex-valued

see the those demo MATLAB scripts for more details.

author : Benyuan Liu
email  : liubenyuan ## gmail ** com
collaborate with : Zhilin Zhang (zhilinzhang ## ieee ** org)

codes:

BSBL_FM.m            ---->    the main algorithm, also called MBSBL-FM in MMV model
BSBL_BO.m            ---->    Zhilin's BSBL-BO algorithm.
demo_smv_real.m      ---->    the real-SMV-block sparse demo
demo_smv_complex.m   ---->    the complex-SMV-block sparse demo
demo_mmv.m           ---->    the real-MMV-block sparse demo
demo_fecg.m          ---->    the demo code for FECG data recovery

data:

demo.mat         ---->    the data for SMV case, contains re, im vectors
signal_01.mat    ---->    FECG datasets used in BSBL by Zhilin
Phi.mat          ---->    the sensing matrix for CS FECG data

