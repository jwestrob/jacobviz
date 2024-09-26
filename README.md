
# Interactive Contact Map Visualization

## Overview

Interactive Contact Map Visualization is a web-based tool that allows users to visualize and analyze protein contact maps alongside 3D protein structures. The application provides interactive heatmaps, sequence alignments, and 3D molecular visualizations to facilitate the exploration of inter-subunit contacts and residue interactions.

## Features

- **Interactive Heatmap:** Visualize contact maps with dynamic highlighting of residues and inter-subunit contacts.
- **3D Protein Structure:** Explore the protein structure in 3D with distinct colors for each subunit.
- **Residue Labels:** Display labels with residue letters and indices next to highlighted residues.
- **Sequence Alignment:** Align input sequences with PDB sequences using the Needleman-Wunsch algorithm for accurate residue mapping.
- **Dynamic Controls:** Generate contact maps from CIF files with customizable distance thresholds and chain selections.

## Installation

### Prerequisites

- **Node.js** (v14 or higher)
- **npm** (Node Package Manager)

### Steps

1. **Clone the Repository**

   ```bash
   git clone https://github.com/yourusername/contact-map-visualization.git
   cd contact-map-visualization
   ```

2. **Install Dependencies**

   ```bash
   npm install
   ```

3. **Start the Application**

   ```bash
   npm start
   ```

4. **Access the Application**

   Open your browser and navigate to `http://localhost:3000`.

## Usage

1. **Select Contact Map**

   - Use the dropdown menu to select a contact map JSON file.
   - The sequences and interactive heatmap will load automatically.

2. **Interact with Heatmap**

   - Hover over cells in the heatmap to highlight corresponding residues in the sequences and 3D structure.
   - Click on residues in the sequence to manually highlight them.

3. **3D Protein Visualization**

   - Explore the protein structure in 3D.
   - Highlighted residues will display labels with residue letters and indices.

4. **Generate Contact Map from CIF File**

   - Select chains and set a distance threshold.
   - Click "Generate Contact Map" to visualize contacts based on the CIF file.

## Project Structure

- `index.html` - Main HTML file.
- `styles.css` - CSS styles for the application.
- `script.js` - Main JavaScript file handling UI interactions and visualization.
- `server.js` - Main JavaScript file handling server and web worker interactions.
- `pdbWorker.js` - Web Worker script for parsing PDB files and sequence alignment.
- `data/` - Directory containing contact map JSON files.

## Technologies Used

- **D3.js:** For creating dynamic and interactive data visualizations.
- **3Dmol.js:** For rendering 3D molecular structures.
- **Web Workers:** For offloading heavy computations like sequence alignment to improve performance.
- **Node.js & Express:** Backend server to serve data files and handle API requests.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

For any inquiries or support, please contact Jacob West-Roberts [jacob@tatta.bio](mailto:jacob@tatta.bio).
