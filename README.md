
	 ______   ______  __   ______   ______   _____    ______    
	/\  ___\ /\  == \/\ \ /\  ___\ /\  __ \ /\  __-. /\  ___\   
	\ \  __\ \ \  _-/\ \ \\ \ \____\ \ \/\ \\ \ \/\ \\ \  __\   
	 \ \_____\\ \_\   \ \_\\ \_____\\ \_____\\ \____- \ \_____\ 
	  \/_____/ \/_/    \/_/ \/_____/ \/_____/ \/____/  \/_____/ 
	                                                            

_MATLAB code to analyze projects of the research group:_
 _"Cellular excitability and neuronal network dynamics", ICM, Paris_

## Goals

 * Increase robustness and automatization of our analyses.
 * Develop and share code and analysis pipelines that are common between research projects.
 * Make analysis pipelines and parameters visible and comparable between projects.
 * Reduce repetition/duplication of code (and work)
 * Provide backup and version control, and a way to update between computers/servers.
 * Provide a single solution for interfacing with in-house developed software (MUSE) and externally developed software ([Spyking-Circus](https://github.com/spyking-circus)).
 * Facilitate documentation and publication of methodology.

## Design

 * A single [MATLAB struct](https://www.mathworks.com/help/matlab/ref/struct.html) contains all settings.
 * This settings struct is passed to all functions, so that the analysis script stays clean.
 * Data, results and figures are saved in project-specific directories
 * Important steps in the analyses can be either read if done earlier, or run again.
 * Processed data follows as much as possible [FieldTrip](https://github.com/fieldtrip/fieldtrip) conventions, making it easier to use FieldTrip at any point of the analysis pipeline. In fact, much of the code uses functions from [FieldTrip](https://github.com/fieldtrip/fieldtrip), and EpiCode functions can be considered higher-level wrapper functions for [FieldTrip](https://github.com/fieldtrip/fieldtrip) functions. 
 * Because many of our analyses revolve around events identified manually with MUSE, another [MATLAB struct](https://www.mathworks.com/help/matlab/ref/struct.html) is created, which can be used to keep track of event timings, as well as interface with the MUSE software.

## Organization

EpiCode is organized as follows:

 * [shared](shared): Main code, and their documentation, shared between projects 
 * [external](external): External code, i.e. not made by me, and not within the EpiCode GNU license
 * [development](development): Code in development, and old code that is saved, just in case
 * [projects](projects): Contains the analysis projects. These are organized as follows:

 Every analysis project directory contains:
 * A single settings file that contains all the parameters for all the analyses
 * An analysis pipeline, i.e. an analysis script that calls EpiCode (and other) fuctions, and i.e. loops over patient.
 * An [R script](https://www.r-project.org/) for statistical analyses
 * A MATLAB and a bash/SLURM script, for running analyses on the computing cluster

## Dependencies

 * The EpiCode relies heavily on [FieldTrip](https://github.com/fieldtrip/fieldtrip), which should be on your MATLAB path, and regularly updated.
 * Spike analyses relie on [Spyking-Circus](https://github.com/spyking-circus).
 * The creation of markers/annotations in the data is done with in-house developed software (MUSE), and the functions and analysis pipeline reflects this.

## Limitations and considerations

 * Use at your own risk.
 * The code is under heavy development, so but will hopefully stabilize over time.
 * The code is aimed specifically at analysing intercranial data that is recorded from epileptic patients over the course of presurgical evaluation. 
 * No data is or will evert be shared in this repository. 

## Reporting issues

 * Please use the [issues page](https://github.com/stephenwhitmarsh/EpiCode/issues) to report and keep track of bugs and requests.

## License

 * GNU GPL v3, see LICENSE file for more details. [External](external) are shared according to their own license files (I lost some, but will include them ASAP) 
