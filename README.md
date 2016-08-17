
https://zenodo.org/badge/23914/KIT-HYD/SeriesDistance.svg

Documentation of the Series Distance (SD) method

**Contact**: uwe.ehret@kit.edu | simon.seibert@kit.edu

**last modification**: 17.08.2016

This is a quick start documentation of the matlab code on the Series Distance (SD) method described by Ehret and Zehe (2011) and Seibert et al. (2016) to quantify timing and magnitude errors between simulated and observed stream flow. Please note that this version of SD is a **development release** which might be unstable and behave unexpectedly. The SD method comes with ABSOLUTELY NO WARRANTY. You are welcome to modify and redistribute it within the license agreement. The SD method is published under the CreativeComons "BY-NC-SA 4.0" license together with a ready-to-use sample data set. To view a full version of the license agreement please visit https://creativecommons.org/licenses/by-nc-sa/4.0/. 

**Concept:**
The SD method is available in a 'continuous' and an 'event mode' version. 
The continuous version requires a simulated (sim) and a related observed (obs) time series as input and produces a single, 2-d distribution of timing and magnitude errors as output.
The event mode version requires a simulated and a related observed time series plus information on the starts and ends of hydrological 'events' (times of high flow) within the time series. It produces three separate 2-d error distributions of timing and magnitude error for i) rising and ii) falling segments within events, and iii) for the remaining periods of low-flow.

**SD_Analysis_Continuous**
* Input (ascii format)
  - obs: [n,1] matrix with equidistant and NaN-free time series of observed discharge data.
  - sim: [n,1] matrix with equidistant and NaN-free time series of simulated discharge data. 
  - timeseries_splits: [n,1] optional matrix with points in time where the 'obs' and 'sim' time series will be split to increase computational speed.

* Parameters are explained and specified in the parameter block of the program

* Output (all stored within a single binary file .mat)
  - obs:               observed discharge (smoothed)
  - sim:               simulated discharge
  - segs_obs_opt:      optimized segments in obs
  - segs_sim_opt:      optimized segments in sim
  - connectors:        SD connectors
  - e_sd_q_all:        SD magnitude errors
  - e_sd_t_all:        SD timing errors
  - parameters:        all parameters, i.e. error_model, objective function weights, smoothing parameters,  
                       are included in the outputfile

**SD_Analysis_Event**
* Input (ascii format)
  - obs: [n,1] matrix with equidistant and NaN-free time series of observed discharge data.
  - sim: [n,1] matrix with equidistant and NaN-free time series of simulated discharge data. 
  - obs_events: [m,2] matrix with start and end times of events in 'obs'. [m,1]=start, [m,2]=end.
  - sim_events: [m,2] matrix with start and end times of events in 'sim'. [m,1]=start, [m,2]=end. 
  - obs_sim_pairing: [m,2] matrix with events in 'obs_events' and 'sim_events' that belong together.

* Parameters are explained and specified in the parameter block of the program

* Output (all stored within a single binary file .mat)
  - obs:               observed discharge (smoothed)
  - sim:               simulated discharge
  - segs_obs_opt_all:  optimized segments in obs
  - segs_sim_opt_all:  optimized segments in sim
  - connectors:        SD connectors
  - e_sd_q_rise:       SD magnitude errors for rising limbs
  - e_sd_t_rise:       SD timing errors for rising limbs
  - e_sd_q_fall:       SD magnitude errors for falling limbs
  - e_sd_t_fall:       SD timing errors for falling limbs
  - e_sd_lowFlow:      1D errors for low-flow
  - parameters:        all parameters, i.e. error_model, objective function weights, smoothing parameters,  
                       are included in the outputfile
  - error_model:       'relative' or 'standard'
  - seg_raw_statistics: statistics of the original time series e.g. number of peaks, throughs, duration of rising/ falling segments
  - seg_opt_statistics: statistics of the coarse-grained time series e.g. number of peaks, throughs, duration of rising/ falling segments

**Getting started:** A test data set is provided in folder \data. To run either SD_Analysis_Continuous or SD_Analysis_Event, adjust the path settings in the script (section 'specify paths and parameters, read inputs'). If you want to use your own data, copy them into \data and adjust the file name settings in the same section. Set parameters in the same section.

**References:**

Ehret, U., Zehe, E., 2011. Series distance - An intuitive metric to quantify hydrograph similarity in terms of occurrence, amplitude and timing of hydrological events. Hydrol. Earth Syst. Sci. 15, 877–896. doi:10.5194/hess-15-877-2011

Seibert, S. P., Ehret, U., and Zehe, E., 2016: Disentangling timing and amplitude errors in streamflow simulations, Hydrol. Earth Syst. Sci. XXXX
