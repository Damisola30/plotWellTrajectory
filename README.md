# PlotWellTrajectory

## Overview

**PlotWellTrajectory** is a MATLAB-based application designed specifically for drilling engineers to visualize and analyze well trajectories. The tool streamlines the process of plotting well paths, making it easier to interpret and communicate drilling plans and survey data.

## Features

- **Intuitive User Interface**: Simplifies the process of loading, visualizing, and analyzing well trajectory data.
- **2D and 3D Visualization**: Supports both planar and spatial representations of well paths.
- **Survey Data Input**: Accepts standard directional survey data formats (measured depth, inclination, azimuth).
- **Flexible Plotting Options**: Allows users to plot vertical section, plan view, and 3D perspective.
- **Engineering-Focused**: Built with the needs of drilling engineers, facilitating rapid interpretation and reporting.

## Prerequisites

- **MATLAB**: This project requires MATLAB to be installed on your system. 
  - Recommended: MATLAB R2018b or newer.
  - Some features may require the MATLAB Plotting Toolbox.
- **Directional Survey Data**: Your well trajectory data should be available in a format compatible with the application (CSV, Excel, or MATLAB table with columns for MD, Inclination, Azimuth).

## Getting Started

1. **Clone or Download the Repository**
    ```bash
    git clone https://github.com/Damisola30/plotWellTrajectory.git
    ```
    Or download the ZIP and extract it.

2. **Open in MATLAB**
    - Launch MATLAB.
    - Add the project folder to your MATLAB path:
      ```matlab
      addpath('path_to_cloned_folder/plotWellTrajectory')
      ```

3. **Run the Application**
    - Navigate to the folder in MATLAB.
    - Run the main script or function, for example:
      ```matlab
      PlotWellTrajectory
      ```
    - Follow the prompts to load your survey data and generate plots.

## Usage

- **Loading Data**: The application will prompt you to select a data file or paste survey values directly. Supported formats include CSV and Excel.
- **Plotting Options**: Choose between 2D (vertical section, plan view) and 3D plots for comprehensive visualization.
- **Customization**: Adjust plot settings such as color, marker size, and axis limits as needed.
- **Exporting Plots**: Save figures in standard image formats for reports and presentations.

## Example

```matlab
% Example usage for a typical workflow

% Load data (replace with your file path)
data = readtable('well_survey.csv');

% Plot 3D trajectory
PlotWellTrajectory(data.MD, data.Inclination, data.Azimuth);

% Customize your plot (see documentation within the app for more options)
```

## Documentation

For detailed documentation, refer to the provided PDF: [`PlotWellTrajectory_DocumentationGroup04.pdf`](PlotWellTrajectory_DocumentationGroup04.pdf)

## Support

For questions, issues, or suggestions, please open an issue on the [GitHub repository](https://github.com/Damisola30/plotWellTrajectory/issues).

---

**Note**: This tool is intended for educational and engineering use. Always verify results and consult with your engineering team before making operational decisions.
````
