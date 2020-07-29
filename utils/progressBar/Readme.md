# Parfor progress monitor
<img align="right" src="https://github.com/DylanMuir/ParforProgMon/raw/master/progress_bar.png" />

## A Java-based `Matlab` class for progress monitoring during a `parfor` loop

### Usage
Begin by creating a parallel pool.
 
Then construct a ParforProgMon object:

    ppm = ParforProgMon(strWindowTitle, nNumIterations <, nProgressStepSize, nWidth, nHeight>);
 
 `strWindowTitle` is a string containing the title of the progress bar
  window. `nNumIterations` is an integer with the total number of
  iterations in the loop.
 
#### Optional arguments
  `nProgressStepSize` specifies to update the progress bar every time this
  number of steps passes. `nWidth` and `nHeight` specify the size of the
  progress window.
 
Within the `parfor` loop:

    parfor (nIndex = 1:nNumIterations)
       ppm.increment();
    end

### Credits
[Parfor Progress monitor](https://www.mathworks.com/matlabcentral/fileexchange/24594-parfor-progress-monitor)

[Parfor Progress monitor v2](https://www.mathworks.com/matlabcentral/fileexchange/31673-parfor-progress-monitor-v2)
