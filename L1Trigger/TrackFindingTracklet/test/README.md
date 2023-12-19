## Overview

In this folder of this branch, there are numerous plotting scripts that are useful for debugging different versions of the Kalman Filter. Below is a table that names the code, gives and example, and a short description. You can find each code in this folder. Each code has some instructions on how to use it in the beginning.

## Codes

| Code Name | Example | Description |
| ----------| ------- | ----------- |
| res_vs_eta_seeds_and_num.C |  ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/a18cd609-6d7d-4d0c-b757-346ddb8fcbbc) ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/b7045922-0152-4ce8-b356-a7d4d6e15314) | Creates two plots. One of them is resolution of z0 plotted against absolute eta separated by seed type. The second is the number of tracks against eta per seed type. |
| res_vs_eta_nstub_and_num.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/20e0927b-4db8-4080-ad0a-85fd03634c9c) ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/4a0aaffb-17d0-4b82-a078-d5921f782efe) | Creates two plots. One of them is resolution of z0 plotted against absolute eta separated by number of stubs. The second is the number of tracks against eta per stub number. |
| overlayHists_resVsEta_3plots.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/373871b7-3352-406e-9e33-afae8a6c2d1a) | Makes 4 resVsEta plots for the four track parameters. Each is named accordingly automatically after the custom name you give it. You can also set each max value manually in an array.
| nstubs_vs_eta.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/334a7a9b-5b07-4b77-8944-862d28dce82f) | Plots Average or RMS of # of stubs against absolute eta|
| ZRplot.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/119d04ba-fbc7-4b00-a234-b99f9d8b591d) | This plot is useful for seeing a visual of how a matched track could have issues in relation to a tracking particle and teh stubs. If (like in this plot) the matched track is showing obviou issues like not passing through stubs in multiple layers, then that's a great indicator of where to look next to fix your bug. |
| ZRplot_badEta_finder.C | ![image](https://github.com/cms-L1TK/cmssw/assets/71595540/50fa9aa2-5d1c-43b6-a3fe-5571cfc7ea85) | This code is a supplement to ZRplot.C by showing which events and tp indicies have particularly bad eta residuals. You can set thresholds, and it will list all of those which fall within a certain eta residual range. Then you can take those and plug them into ZRplot.C to see what the innerworkings of the tracking look like in situations when the residual it particularly bad or good. |
