# flow_evolution_and_solute_dispersion_in_the_swash_of_regular_waves
Code and data used in the manuscript "Flow evolution and solute dispersion in the swash of regular waves" submitted to JGR: Oceans.

1_camera_preprocessing/

  * A00_Smooth_background.tiff - Beach background image with still water.
  * A00_SWL.mat - Location of still water line.
  * Calibration_Results_072424.mat - Camera calibration results.
  * camera_image_preprocess.m - Code for image rectification/processing.
  * cc_smooth.mat - Calibration coefficient to convert between pixels and mm.

2_particle_tracking/
  * ea_cycles.m - Compute individual wave cycles.
  * find_particles_1.m - STEP 1: Find particle locations.
  * build_tracks_2.m - STEP 2: Connect particle locations into tracks.
  * process_track_details_3.m - STEP 3: Compute and save track details.
  * save_sptv_4.m - STEP 4: Save SPTV data.
  * ens_av_ptv_scatter.m - Scatter and save SPTV data.
  * ensemble_avg_frames.m - Use shoreline and still water line to cycle align frames.
  * track_details_v2_particle_centroids.m - Function to compute track details.
  * track_parts.m - Function to connect particles into tracks.
  * velocity_fn.m - Velocity convolution function.

3_shoreline_tracking/
  * ea_cycles.m - Compute individual wave cycles.
  * ensemble_avg_frames.m - Use shoreline and still water line to cycle align frames.
  * optical_shoreline_final.m - Compute optical shoreline location.
  * combine_ea_shoreline.m - Compute enhanced shoreline with SPTV results from two wave trials.

4_sensors/
  * combine_trials.m - Combine two SPTV trials.
  * ensemble_average_sensors_depth.m - Ensemble average and combine the depth data from wave gauges across two independent trials.
  * ensemble_average_sensors_vel.m - Ensemble average and combine the velocity data from velocimeters across two independent trials.
  * ea_cycles.m - Compute individual wave cycles,
  * ensemble_average_sensors.m - Ensemble average depth and velocity data from sensors and SPTV point measurements,
  * ensemble_avg_frames.m - Use shoreline and still water line to cycle align frames,
  * WG_locs.mat - Location of wave gauges in camera FOV
  * Raw sensor data is too large for gitHub, but can be downloaded from doi:10.17603/ds2-vy34-ad78
      * H = 0.1m: "SMOOTH BED/Angle00/Case03R/RAW/Sensors/RAW_SMOOTH_ANGLE00_Case01R.mat"
      * H = 0125m: "SMOOTH BED/Angle00/Case03R/RAW/Sensors/RAW_SMOOTH_ANGLE00_Case02R.mat"
      * H = 0.15m: "SMOOTH BED/Angle00/Case03R/RAW/Sensors/RAW_SMOOTH_ANGLE00_Case03R.mat"

6_depth_field/
  * depth_field.m - Interpolate depth field in the swash zone from local sensors,
  * WG_locs.mat - Location of wave gauges in camera FOV

7_dye/
  * A00_Smooth1 - First dye drop data,
  * A00_Smooth2 - Second dye drop data,
  * A00_Smooth3 - Thrid dye drop data,
  * A00_Smooth4 - Fourth dye drop data,
  * dye_cloud_1.tif - Sample dye cloud image,
  * shoreline_dye.mat - Dye trial shoreline data

8_model/
  * Antuono2010.m - Code to solve Antuono2010 swash zone model

9_Figures/
  * Figure_2.m
  * Figure_4.m
  * Figures_5_6_7.m
  * Fiugre_8.m
  * Figures_9_10_11.m
  * Figure_12_13.m
  * Figures_14_15_16/
    * dye_cloud_4bcde.m - Plot Figure 4b-e,
    * dye_compensated.m - Plot Figure S8,
    * dye_figs_final.m - Plot Figure 4a,
    * dye_image.m - Plot sample dye image,
    * ea_cycles.m - Compute individual wave cycles,
    * ensemble_avg_frames_dye.m - Use shoreline and still water line to cycle align frames,

Results/
  * background.mat - Camera background image (average of 100 random frames during particle trial),
  * combined_shoreline_ea.mat - Shoreline computed from both optical shoreline and enhanced by SPTV results, H = 0.15m
  * connected_particles_byframe.mat - All particles identified in each frame with particle identifies noting connected particles, H = 0.15m
  * connected_particles.mat - Each connected particle track, H = 0.15m
  * depth_ensembel_average.mat - Ensemble average depths recorded from wave gagues, H = 0.15m
  * ea_detph_field.mat - Ensemble average depth field in the swash zone, H = 0.15m
  * EA_uxt.mat - Particle positions and velocities at cycle aligned times [0-2]s, H = 0.15m
  * particle_locations.mat - All identified particle locations, H = 0.15m
  * ptv_xyuv.mat - SPTV data: each frame lists identified velocities at specific positions, H = 0.15m
  * shoreline.mat - Optical shoreline results, H = 0.15m
  * time_trim.mat - Laboratory time of each camera frame (in seconds), H = 0.15m
  * track_details_px_f.mat - All track details in [pixels and frames], H = 0.15m
  * velocity_ensemble_average.mat - Ensemble average velocities at each sensor/POI, H = 0.15m
  * Additional data for other wave cases available from  doi:10.17603/ds2-vy34-ad78

Sample_Images/
  * Several sample SPTV images (rectified)
