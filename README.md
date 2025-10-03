# flow_evolution_and_solute_dispersion_in_the_swash_of_regular_waves
Code and data used in the manuscript "Flow evolution and solute dispersion in the swash of regular waves" submitted to GRL.

1_camera_preprocessing/

  A00_Smooth_background.tiff - Beach background image with still water,
  A00_SWL.mat - Location of still water line,
  Calibration_Results_072424.mat - Camera calibration results,
  camera_image_preprocess.m - Code for image rectification/processing,
  cc_smooth.mat - Calibration coefficient to convert between pixels and mm

2_particle_tracking/
  ea_cycles.m - Compute individual wave cycles,
  find_particles_1.m - STEP 1: Find particle locations,
  build_tracks_2.m - STEP 2: Connect particle locations into tracks,
  process_track_details_3.m - STEP 3: Compute and save track details,
  save_sptv_4.m - STEP 4: Save SPTV data,
  ens_av_ptv_scatter.m - Scatter and save SPTV data,
  ensemble_avg_frames.m - Use shoreline and still water line to cycle align frames,
  track_details_v2_particle_centroids.m - Function to compute track details,
  track_parts.m - Function to connect particles into tracks,
  velocity_fn.m - Velocity convolution function

3_shoreline_tracking/
  ea_cycles.m - Compute individual wave cycles,
  ensemble_avg_frames.m - Use shoreline and still water line to cycle align frames,
  optical_shoreline_final.m - Compute optical shoreline location,
  new_shoreline_v2.m - Compute enhanced shoreline with SPTV results

4_sensors/
  datalogger.mat - All data saved to datalogger (wave gauges and paddle trigger),
  ea_cycles.m - Compute individual wave cycles,
  ensemble_average_sensors.m - Ensemble average depth and velocity data from sensors and SPTV point measurements,
  ensemble_avg_frames.m - Use shoreline and still water line to cycle align frames,
  WG_locs.mat - Location of wave gauges in camera FOV
  Velocimeters_Smooth_A00_C03_T01.mat - ADV data.  File is too large for gitHub but can be downloaded from doi:10.17603/DS2-VY34-AD7

6_depth_field/
  depth_field.m - Interpolate depth field in the swash zone from local sensors,
  WG_locs.mat - Location of wave gauges in camera FOV

7_dye/
  A00_Smooth1 - First dye drop data,
  A00_Smooth2 - Second dye drop data,
  A00_Smooth3 - Thrid dye drop data,
  A00_Smooth4 - Fourth dye drop data,
  dye_cloud_1.tif - Sample dye cloud image,
  shoreline_dye.mat - Dye trial shoreline data

8_model/
  Antuono2010.m - Code to solve Antuono2010 swash zone model

9_paper_figures/
  Fig_Dye/
    dye_cloud_4bcde.m - Plot Figure 4b-e,
    dye_compensated.m - Plot Figure S8,
    dye_figs_final.m - Plot Figure 4a,
    dye_image.m - Plot sample dye image,
    ea_cycles.m - Compute individual wave cycles,
    ensemble_avg_frames_dye.m - Use shoreline and still water line to cycle align frames,
  Fig_ref/
    Figure1.m - Plot Figure 1,
    Figure_S5.m - Plot Figure S5,
  FigN_Contours/
    contours_data.m - Figure 3 data panels,
    contours_model.m - Figure 3 model panels,
    crameri_colormaps - Fabio Crameri's colormaps for plotting,
    slanCM_Data.mat - Colormaps from https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap,
    slanCM.m - Colormap function from https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap ,
    thiel_sen.m - Thiel-Sen estimator function,
  FigN_Points/
    figure_point_measurements.m - Plot Figure 2,
    plot_alpha_ref.m - Function for plotting Figure 2 alpha,
    plot_depth_ref.m - Function for plotting Figure 2 depth,
    plot_velocity_ref.m - Function for plotting Figure 2 velocity,
    sw.mat - Wave gauge still water depths,
    WG_locs.mat - Wave gauge locations in camera FOV

Results/
  background.mat - Camera background image (average of 100 random frames during particle trial),
  combined_shoreline_ea.mat - Shoreline computed from both optical shoreline and enhanced by SPTV results,
  connected_particles_byframe.mat - All particles identified in each frame with particle identifies noting connected particles,
  connected_particles.mat - Each connected particle track,
  depth_ensembel_average.mat - Ensemble average depths recorded from wave gagues,
  ea_detph_field.mat - Ensemble average depth field in the swash zone,
  EA_uxt.mat - Particle positions and velocities at cycle aligned times [0-2]s,
  particle_locations.mat - All identified particle locations,
  ptv_xyuv.mat - SPTV data: each frame lists identified velocities at specific positions,
  shoreline.mat - Optical shoreline results,
  time_trim.mat - Laboratory time of each camera frame (in seconds),
  track_details_px_f.mat - All track details in [pixels and frames],
  velocity_ensemble_average.mat - Ensemble average velocities at each sensor/POI.

Sample_Images/
  Several sample SPTV images (rectified)
