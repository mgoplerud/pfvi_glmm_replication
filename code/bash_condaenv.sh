# Create the conda environment
conda create -n "python_modern" python==3.9.1
# Install a few packages needed
source activate python_modern
pip install polyagamma bagpipes numpy scipy pandas --upgrade

