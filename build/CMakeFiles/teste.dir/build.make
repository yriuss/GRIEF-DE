# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/adriel/teste/GRIEF-DE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/adriel/teste/GRIEF-DE/build

# Include any dependencies generated for this target.
include CMakeFiles/teste.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/teste.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/teste.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/teste.dir/flags.make

CMakeFiles/teste.dir/main.cu.o: CMakeFiles/teste.dir/flags.make
CMakeFiles/teste.dir/main.cu.o: ../main.cu
CMakeFiles/teste.dir/main.cu.o: CMakeFiles/teste.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adriel/teste/GRIEF-DE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/teste.dir/main.cu.o"
	/usr/local/cuda-11.2/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/teste.dir/main.cu.o -MF CMakeFiles/teste.dir/main.cu.o.d -x cu -dc /home/adriel/teste/GRIEF-DE/main.cu -o CMakeFiles/teste.dir/main.cu.o

CMakeFiles/teste.dir/main.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/teste.dir/main.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/teste.dir/main.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/teste.dir/main.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target teste
teste_OBJECTS = \
"CMakeFiles/teste.dir/main.cu.o"

# External object files for target teste
teste_EXTERNAL_OBJECTS =

CMakeFiles/teste.dir/cmake_device_link.o: CMakeFiles/teste.dir/main.cu.o
CMakeFiles/teste.dir/cmake_device_link.o: CMakeFiles/teste.dir/build.make
CMakeFiles/teste.dir/cmake_device_link.o: ../lib/libGRIEF_proj.so
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_gapi.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_stitching.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_alphamat.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_aruco.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_bgsegm.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_bioinspired.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_ccalib.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudabgsegm.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudafeatures2d.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudaobjdetect.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudastereo.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_dnn_objdetect.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_dnn_superres.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_dpm.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_face.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_freetype.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_fuzzy.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_hdf.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_hfs.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_img_hash.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_intensity_transform.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_line_descriptor.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_mcc.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_quality.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_rapid.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_reg.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_rgbd.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_saliency.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_sfm.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_stereo.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_structured_light.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_phase_unwrapping.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_superres.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_surface_matching.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_tracking.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_highgui.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_datasets.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_plot.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_text.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_videostab.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_videoio.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudaoptflow.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudalegacy.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudawarping.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_optflow.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_wechat_qrcode.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_xfeatures2d.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_ml.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_shape.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_ximgproc.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_video.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_dnn.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_xobjdetect.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_imgcodecs.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_objdetect.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_calib3d.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_features2d.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_flann.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_xphoto.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_photo.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudaimgproc.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudafilters.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_imgproc.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudaarithm.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_core.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/local/lib/libopencv_cudev.so.4.5.2
CMakeFiles/teste.dir/cmake_device_link.o: /usr/lib/x86_64-linux-gnu/libpython3.8.so
CMakeFiles/teste.dir/cmake_device_link.o: CMakeFiles/teste.dir/dlink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/adriel/teste/GRIEF-DE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA device code CMakeFiles/teste.dir/cmake_device_link.o"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/teste.dir/dlink.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/teste.dir/build: CMakeFiles/teste.dir/cmake_device_link.o
.PHONY : CMakeFiles/teste.dir/build

# Object files for target teste
teste_OBJECTS = \
"CMakeFiles/teste.dir/main.cu.o"

# External object files for target teste
teste_EXTERNAL_OBJECTS =

teste: CMakeFiles/teste.dir/main.cu.o
teste: CMakeFiles/teste.dir/build.make
teste: ../lib/libGRIEF_proj.so
teste: /usr/local/lib/libopencv_gapi.so.4.5.2
teste: /usr/local/lib/libopencv_stitching.so.4.5.2
teste: /usr/local/lib/libopencv_alphamat.so.4.5.2
teste: /usr/local/lib/libopencv_aruco.so.4.5.2
teste: /usr/local/lib/libopencv_bgsegm.so.4.5.2
teste: /usr/local/lib/libopencv_bioinspired.so.4.5.2
teste: /usr/local/lib/libopencv_ccalib.so.4.5.2
teste: /usr/local/lib/libopencv_cudabgsegm.so.4.5.2
teste: /usr/local/lib/libopencv_cudafeatures2d.so.4.5.2
teste: /usr/local/lib/libopencv_cudaobjdetect.so.4.5.2
teste: /usr/local/lib/libopencv_cudastereo.so.4.5.2
teste: /usr/local/lib/libopencv_dnn_objdetect.so.4.5.2
teste: /usr/local/lib/libopencv_dnn_superres.so.4.5.2
teste: /usr/local/lib/libopencv_dpm.so.4.5.2
teste: /usr/local/lib/libopencv_face.so.4.5.2
teste: /usr/local/lib/libopencv_freetype.so.4.5.2
teste: /usr/local/lib/libopencv_fuzzy.so.4.5.2
teste: /usr/local/lib/libopencv_hdf.so.4.5.2
teste: /usr/local/lib/libopencv_hfs.so.4.5.2
teste: /usr/local/lib/libopencv_img_hash.so.4.5.2
teste: /usr/local/lib/libopencv_intensity_transform.so.4.5.2
teste: /usr/local/lib/libopencv_line_descriptor.so.4.5.2
teste: /usr/local/lib/libopencv_mcc.so.4.5.2
teste: /usr/local/lib/libopencv_quality.so.4.5.2
teste: /usr/local/lib/libopencv_rapid.so.4.5.2
teste: /usr/local/lib/libopencv_reg.so.4.5.2
teste: /usr/local/lib/libopencv_rgbd.so.4.5.2
teste: /usr/local/lib/libopencv_saliency.so.4.5.2
teste: /usr/local/lib/libopencv_sfm.so.4.5.2
teste: /usr/local/lib/libopencv_stereo.so.4.5.2
teste: /usr/local/lib/libopencv_structured_light.so.4.5.2
teste: /usr/local/lib/libopencv_phase_unwrapping.so.4.5.2
teste: /usr/local/lib/libopencv_superres.so.4.5.2
teste: /usr/local/lib/libopencv_surface_matching.so.4.5.2
teste: /usr/local/lib/libopencv_tracking.so.4.5.2
teste: /usr/local/lib/libopencv_highgui.so.4.5.2
teste: /usr/local/lib/libopencv_datasets.so.4.5.2
teste: /usr/local/lib/libopencv_plot.so.4.5.2
teste: /usr/local/lib/libopencv_text.so.4.5.2
teste: /usr/local/lib/libopencv_videostab.so.4.5.2
teste: /usr/local/lib/libopencv_videoio.so.4.5.2
teste: /usr/local/lib/libopencv_cudaoptflow.so.4.5.2
teste: /usr/local/lib/libopencv_cudalegacy.so.4.5.2
teste: /usr/local/lib/libopencv_cudawarping.so.4.5.2
teste: /usr/local/lib/libopencv_optflow.so.4.5.2
teste: /usr/local/lib/libopencv_wechat_qrcode.so.4.5.2
teste: /usr/local/lib/libopencv_xfeatures2d.so.4.5.2
teste: /usr/local/lib/libopencv_ml.so.4.5.2
teste: /usr/local/lib/libopencv_shape.so.4.5.2
teste: /usr/local/lib/libopencv_ximgproc.so.4.5.2
teste: /usr/local/lib/libopencv_video.so.4.5.2
teste: /usr/local/lib/libopencv_dnn.so.4.5.2
teste: /usr/local/lib/libopencv_xobjdetect.so.4.5.2
teste: /usr/local/lib/libopencv_imgcodecs.so.4.5.2
teste: /usr/local/lib/libopencv_objdetect.so.4.5.2
teste: /usr/local/lib/libopencv_calib3d.so.4.5.2
teste: /usr/local/lib/libopencv_features2d.so.4.5.2
teste: /usr/local/lib/libopencv_flann.so.4.5.2
teste: /usr/local/lib/libopencv_xphoto.so.4.5.2
teste: /usr/local/lib/libopencv_photo.so.4.5.2
teste: /usr/local/lib/libopencv_cudaimgproc.so.4.5.2
teste: /usr/local/lib/libopencv_cudafilters.so.4.5.2
teste: /usr/local/lib/libopencv_imgproc.so.4.5.2
teste: /usr/local/lib/libopencv_cudaarithm.so.4.5.2
teste: /usr/local/lib/libopencv_core.so.4.5.2
teste: /usr/local/lib/libopencv_cudev.so.4.5.2
teste: /usr/lib/x86_64-linux-gnu/libpython3.8.so
teste: CMakeFiles/teste.dir/cmake_device_link.o
teste: CMakeFiles/teste.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/adriel/teste/GRIEF-DE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CUDA executable teste"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/teste.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/teste.dir/build: teste
.PHONY : CMakeFiles/teste.dir/build

CMakeFiles/teste.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/teste.dir/cmake_clean.cmake
.PHONY : CMakeFiles/teste.dir/clean

CMakeFiles/teste.dir/depend:
	cd /home/adriel/teste/GRIEF-DE/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/adriel/teste/GRIEF-DE /home/adriel/teste/GRIEF-DE /home/adriel/teste/GRIEF-DE/build /home/adriel/teste/GRIEF-DE/build /home/adriel/teste/GRIEF-DE/build/CMakeFiles/teste.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/teste.dir/depend

