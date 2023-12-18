## Overview

In this folder of this branch, there are numerous plotting scripts that are useful for debugging different versions of the Kalman Filter. Below is a table that names the code, gives and example, and a short description. You can find each code in this folder. Each code has some instructions on how to use it in the beginning.

## Codes

| Code Name | Example | Description |
| ----------| ------- | ----------- |
| res_vs_eta_seeds_and_num.C |  ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/95d764fd-d921-48aa-aeb5-23589bfe8e7a) ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/b7045922-0152-4ce8-b356-a7d4d6e15314) | Creates two plots. One of them is resolution of z0 plotted against absolute eta separated by seed type. The second is the number of tracks against eta per seed type. |
| res_vs_eta_nstub_and_num.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/20e0927b-4db8-4080-ad0a-85fd03634c9c) ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/4a0aaffb-17d0-4b82-a078-d5921f782efe) | Creates two plots. One of them is resolution of z0 plotted against absolute eta separated by number of stubs. The second is the number of tracks against eta per stub number. |
| overlayHists_resVsEta_3plots.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/373871b7-3352-406e-9e33-afae8a6c2d1a) | Makes 4 resVsEta plots for the four track parameters. Each is named accordingly automatically after the custom name you give it. You can also set each max value manually in an array.
| | | |
