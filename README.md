             _            _        _           _             _            _            _      
            /\ \         /\ \     /\ \       /\ \           /\ \         /\ \         /\ \    
           /  \ \       /  \ \    \ \ \     /  \ \         /  \ \       /  \ \____   /  \ \   
          / /\ \ \     / /\ \ \   /\ \_\   / /\ \ \       / /\ \ \     / /\ \_____\ / /\ \ \  
         / / /\ \_\   / / /\ \_\ / /\/_/  / / /\ \ \     / / /\ \ \   / / /\/___  // / /\ \_\ 
        / /_/_ \/_/  / / /_/ / // / /    / / /  \ \_\   / / /  \ \_\ / / /   / / // /_/_ \/_/ 
       / /____/\    / / /__\/ // / /    / / /    \/_/  / / /   / / // / /   / / // /____/\    
      / /\____\/   / / /_____// / /    / / /          / / /   / / // / /   / / // /\____\/    
     / / /______  / / /   ___/ / /__  / / /________  / / /___/ / / \ \ \__/ / // / /______    
    / / /_______\/ / /   /\__\/_/___\/ / /_________\/ / /____\/ /   \ \___\/ // / /_______\   
    \/__________/\/_/    \/_________/\/____________/\/_________/     \/_____/ \/__________/   


Analysis scripts and functions for projects of the research group:
_"Cellular excitability and neuronal network dynamics", ICM, Paris_

## Goals
 * Make analyses easier by providing high-level functions for common analyses on our data.
 * Increase robustness and automatization of our analyses.
 * Develop and share code and analysis pipelines that are common between research projects.
 * Make analysis pipelines and parameters visible and comparable between projects.
 * Reduce repetition/duplication of code (and work)
 * Backup and version control, and a way to update between computers/servers.
 * Single solution for interfacing with in-house (MUSE) and external software ([Spyking-Circus](https://github.com/spyking-circus)).
 * Synchronization between different recording systems and their annotation systems (NeuraLynx, MicroMed and DeltaMed)
 * Facilitate use of [FieldTrip](https://github.com/fieldtrip/fieldtrip).
 * Facilitate documentation and publication of methodology.

## Design

 * A single [MATLAB struct](https://www.mathworks.com/help/matlab/ref/struct.html) contains all settings.
 * This settings struct is passed to all functions, so that the analysis script stays clean.
 * Data, results and figures are saved in project-specific directories
 * Important steps in the analyses are read from file if done earlier, or forced to run again.
 * Processed data follows as much as possible [FieldTrip](https://github.com/fieldtrip/fieldtrip) conventions, making it easier to use FieldTrip at any point of the analysis pipeline. In fact, much of the code uses functions from [FieldTrip](https://github.com/fieldtrip/fieldtrip), and EpiCode functions can be considered higher-level wrapper functions for [FieldTrip](https://github.com/fieldtrip/fieldtrip) functions. 
 * Because many of our analyses revolve around events identified manually with MUSE, those are read into a [MATLAB struct](https://www.mathworks.com/help/matlab/ref/struct.html) which is used to keep track of timings, and can be used to read/write markers created/read by MUSE.

## Organization

EpiCode is organized as follows:

 * [shared](shared) contains the main code, and their documentation, shared between projects 
 * [external](external) contains code not made by me, and not within the EpiCode GNU license
 * [development](development) contains code in development, and old code, just in case
 * [projects](projects) contains analysis projects. These are organized as follows.

 Every [analysis project directory](projects) contains:
 * A single settings file that contains all parameters
 * An analysis script that calls EpiCode (and other) MATLAB fuctions, and i.e. loops over patients.
 * An [R script](https://www.r-project.org/) for statistical analyses
 * A MATLAB and bash/SLURM script, for running analyses on the computing cluster

## Dependencies

 * EpiCode relies heavily on [FieldTrip](https://github.com/fieldtrip/fieldtrip), which should be on your MATLAB path, and regularly updated.
 * Spike analyses rely on [Spyking-Circus](https://github.com/spyking-circus).
 * The creation of markers/annotations in the data is done with in-house developed software (MUSE), and the functions and analysis pipeline reflects this.

## Limitations and considerations

 * Use at your own risk.
 * The code is under heavy development, and will be (but will hopefully stabilize over time).
 * The code is aimed specifically at analysing intercranial data that is recorded from epileptic patients over the course of presurgical evaluation. 
 * No data is or will ever be shared through this repository. 

## Reporting issues

 * Please use the [issues page](https://github.com/stephenwhitmarsh/EpiCode/issues) to report and keep track of bugs and requests.

## License

 * GNU GPL v3, see LICENSE file for more details. [External](external) are shared according to their own license files (I lost some, but will include them ASAP) 
