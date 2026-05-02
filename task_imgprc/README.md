## Overview
	Subtitle region extraction from videos using
	Haar 2D-DWT edge detection, with
	multi-frame tracking and saliency-aware preprocessing
	for colored text on bright backgrounds

## Status
	Release

## Contains
	main.c          -- single source: video I/O, DWT pipeline, UI
	raw/            -- input video samples

## Pipeline
	BGR -> {V, salience = V/2 + chroma}  fused single-pass
	salience -> Haar 2D-DWT (1 level)
	per-subband adaptive threshold (LH, HL, HH)
	per-subband morphology (paper figure 10 kernels)
	logical AND of dilated edges
	horizontal closing (line merge)
	contour -> geometry/intensity filter -> NMS
	multi-frame IoU tracker
	per-track screenshot on first stable hit

## Dependencies
	gcc             >= 12
	OpenCV4         (core, imgproc, highgui)
	FFmpeg          (libavformat, libavcodec, libswscale, libavutil)
	GNU Make
	libc, libm

## Build
	make

## Run
	./imgprc raw/AtomicBlonde.mp4
or
	make run

## UI
	Position slider     -- scrub video
	DWT Thresh slider   -- threshold tuning
	Min Confidence      -- confidence gate
	a / h               -- seek -5s
	d / l               -- seek +5s
	A / H               -- seek -30s
	D / L               -- seek +30s
	space / p           -- pause/resume
	c                   -- toggle ROI color enhancement
	ESC                 -- quit

## Output
	screenshots/frame_*.ppm  -- one per new stable subtitle event
