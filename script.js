// script.js

// Ensure the script runs after the page has loaded
document.addEventListener('DOMContentLoaded', function () {
    // API endpoint to get the list of JSON data files
    const dataFilesEndpoint = '/api/data-files';

    // Heatmap dimensions
    let heatmapWidth;
    let heatmapHeight;

    // Initialize a cache object to store fetched JSON data
    const dataCache = {};

    // Initialize variables for protein visualization
    let viewer = null;
    let residueMapping = {};

    // Fetch the list of JSON data files from the server
    fetch(dataFilesEndpoint)
        .then(response => response.json())
        .then(data => {
            const dataFiles = data.data_files;

            // Populate the dropdown menu
            const dataSelect = d3.select('#data-select');
            dataSelect.selectAll('option').remove(); // Clear existing options

            dataSelect.selectAll('option')
                .data(dataFiles)
                .enter()
                .append('option')
                .attr('value', d => d)
                .text(d => d.replace('.json', ''));

            // Set the first option as selected and load it
            if (dataFiles.length > 0) {
                loadData(dataFiles[0]);
                dataSelect.property('value', dataFiles[0]);
            } else {
                dataSelect.append('option')
                    .attr('value', '')
                    .attr('disabled', true)
                    .attr('selected', true)
                    .text('No data files available');
            }

            // Add event listener for dropdown change
            dataSelect.on('change', function(event) {
                const selectedFile = event.target.value;
                loadData(selectedFile);
            });

        })
        .catch(error => {
            console.error('Error fetching data files:', error);
            // Update dropdown to show error
            const dataSelect = d3.select('#data-select');
            dataSelect.selectAll('option').remove();
            dataSelect.append('option')
                .attr('value', '')
                .attr('disabled', true)
                .attr('selected', true)
                .text('Error loading data files');
        });

    // Function to load and render data from a JSON file with caching
    function loadData(jsonFile) {
        // Check if the data is already in the cache
        if (dataCache[jsonFile]) {
            console.log(`Loading ${jsonFile} from cache.`);
            renderVisualization(dataCache[jsonFile]);
            return;
        }

        // If not in cache, fetch the data from the server
        d3.json(`data/${jsonFile}`).then(data => {
            // Store the fetched data in the cache
            dataCache[jsonFile] = data;
            console.log(`Fetched and cached ${jsonFile}.`);
            renderVisualization(data);
        }).catch(error => {
            console.error(`Error loading ${jsonFile}:`, error);
            // Optionally, display an error message in the UI
            d3.select('#heatmap-container')
                .append('div')
                .attr('class', 'error')
                .text(`Error loading ${jsonFile}: ${error.message}`);
        });
    }

    // Function to render the visualization given the data
    function renderVisualization(data) {
        const sequences = parseSequences(data.sequences);
        const sequencePositions = calculateSequencePositions(sequences);
        const matrixData = data.contact_map;
        const matrixSize = matrixData.length;

        // Update heatmap dimensions based on container size
        const containerWidth = document.getElementById('heatmap-container').clientWidth;
        const containerHeight = containerWidth; // Keep it square
        heatmapWidth = containerWidth - 40; // Account for margins
        heatmapHeight = containerHeight - 40;

        // Create scales
        const xScale = d3.scaleBand()
            .domain(d3.range(matrixSize))
            .range([0, heatmapWidth])
            .padding(0.01);

        const yScale = d3.scaleBand()
            .domain(d3.range(matrixSize).reverse()) // Inverted Y-axis
            .range([0, heatmapHeight])
            .padding(0.01);

        const colorScale = d3.scaleSequential(d3.interpolateViridis)
            .domain([d3.min(matrixData.flat()), d3.max(matrixData.flat())]);

        // Render sequences
        renderSequences(sequences, sequencePositions);

        // Render heatmap
        renderHeatmap(matrixData, xScale, yScale, colorScale, sequencePositions);

        // If pdb_id is present, fetch and visualize the protein structure
        if (data.pdb_id) {
            fetchPDBAndRender(data.pdb_id, sequences);
        } else {
            // Clear previous protein visualization if any
            d3.select('#protein-visualizer').html('');
            viewer = null;
            residueMapping = {};
        }
    }

    // Function to parse sequences
    function parseSequences(sequenceInput) {
        return sequenceInput.split('<+>').filter(seq => seq.length > 0).map(seq => seq.trim().toUpperCase());
    }

    // Function to calculate sequence positions within the concatenated sequence
    function calculateSequencePositions(sequences) {
        let sequencePositions = [];
        let currentPos = 0;
        const separatorLength = '<+>'.length;

        sequences.forEach((seq, index) => {
            const start = currentPos;
            const end = currentPos + seq.length - 1;
            sequencePositions.push({
                name: `Sequence ${index + 1}`,
                start: start,
                end: end,
                length: seq.length,
                sequence: seq,
                index: index
            });
            // Update current position (+ separator length)
            currentPos = end + separatorLength + 1;
        });
        return sequencePositions;
    }

    // Function to render sequences
    function renderSequences(sequences, sequencePositions) {
        const sequenceContainer = d3.select('#sequences-container');
        sequenceContainer.html(''); // Clear existing sequences

        // Print input sequences to console for debugging
        console.log("Input Sequences:");
        sequences.forEach((seq, seqIndex) => {
            console.log(`Sequence ${seqIndex}: ${seq}`);
        });

        sequences.forEach((seq, seqIndex) => {
            const seqDiv = sequenceContainer.append('div')
                .attr('id', `sequence-${seqIndex}`)
                .attr('class', 'sequence');

            // Add sequence label
            seqDiv.append('div')
                .attr('class', 'sequence-label')
                .text(`Sequence ${seqIndex + 1}:`)
                .style('margin-bottom', '5px');

            // Create a span for the entire sequence
            const seqSpan = seqDiv.append('div')
                .attr('class', 'sequence-content');

            seq.split('').forEach((residue, resIndex) => {
                seqSpan.append('span')
                    .attr('id', `seq${seqIndex}-res${resIndex}`)
                    .attr('class', 'residue')
                    .text(residue)
                    .on('click', function () {
                        const element = d3.select(this);
                        element.classed('manual-highlight', !element.classed('manual-highlight'));
                    })
                    .on('mouseover', function () {
                        const index = sequencePositions[seqIndex].start + resIndex;
                        highlightResidues(index, index, sequencePositions);
                    })
                    .on('mouseout', function () {
                        clearHighlights();
                    });
            });
        });
    }

    // Function to render heatmap
    function renderHeatmap(matrixData, xScale, yScale, colorScale, sequencePositions) {
        const heatmapContainer = d3.select('#heatmap-container');
        heatmapContainer.html(''); // Clear existing heatmap

        const margin = { top: 20, right: 20, bottom: 20, left: 20 };

        const svg = heatmapContainer.append('svg')
            .attr('width', heatmapWidth + margin.left + margin.right)
            .attr('height', heatmapHeight + margin.top + margin.bottom)
            .call(d3.zoom()
                .scaleExtent([1, 10])
                .on('zoom', zoomed))
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // Define clip path
        svg.append('defs').append('clipPath')
            .attr('id', 'clip')
            .append('rect')
            .attr('width', heatmapWidth)
            .attr('height', heatmapHeight);

        // Create a group for the heatmap cells
        const heatmapGroup = svg.append('g')
            .attr('class', 'heatmap-group')
            .attr('clip-path', 'url(#clip)');

        // Flatten matrix data for easier processing
        const flatData = [];
        for (let y = 0; y < matrixData.length; y++) {
            for (let x = 0; x < matrixData[y].length; x++) {
                flatData.push({
                    value: matrixData[y][x],
                    x: x,
                    y: y
                });
            }
        }

        // Render heatmap cells
        heatmapGroup.selectAll('rect')
            .data(flatData)
            .enter()
            .append('rect')
            .attr('x', d => xScale(d.x))
            .attr('y', d => yScale(d.y))
            .attr('width', xScale.bandwidth())
            .attr('height', yScale.bandwidth())
            .style('fill', d => colorScale(d.value))
            .on('mouseover', function (event, d) {
                handleMouseOver(event, d, sequencePositions, xScale, yScale, svg);
            })
            .on('mouseout', function () {
                handleMouseOut();
            });

        // Zoom function
        function zoomed(event) {
            heatmapGroup.attr('transform', event.transform);
        }

        // For debugging: Print heatmap dimensions
        console.log(`Heatmap Dimensions: Width=${heatmapWidth}, Height=${heatmapHeight}`);
    }

    // Function to handle mouse over event
    function handleMouseOver(event, d, sequencePositions, xScale, yScale, svg) {
        // Highlight corresponding residues
        highlightResidues(d.x, d.y, sequencePositions);

        // Display tooltip at the mouse position
        const tooltip = d3.select('#tooltip');
        tooltip.classed('hidden', false)
            .classed('visible', true)
            .html(`Value: ${d.value.toFixed(4)}<br>X: ${d.x}<br>Y: ${d.y}`)
            .style('left', (event.pageX + 15) + 'px')
            .style('top', (event.pageY - 28) + 'px');
    }

    // Function to handle mouse out event
    function handleMouseOut() {
        clearHighlights();
        d3.select('#tooltip')
            .classed('hidden', true)
            .classed('visible', false);
    }

    // Function to highlight residues
    function highlightResidues(xPos, yPos, sequencePositions) {
        clearHighlights();

        const xMapping = mapIndexToSequence(xPos, sequencePositions);
        const yMapping = mapIndexToSequence(yPos, sequencePositions);

        if (xMapping) {
            d3.select(`#seq${xMapping.seqIndex}-res${xMapping.resIndex}`)
                .classed('hover-highlight', true);
        }
        if (yMapping) {
            d3.select(`#seq${yMapping.seqIndex}-res${yMapping.resIndex}`)
                .classed('hover-highlight', true);
        }

        // Highlight residues in the protein structure
        if (viewer && residueMapping) {
            const residuesToHighlight = [];

            if (xMapping && residueMapping[xPos]) {
                residuesToHighlight.push(...residueMapping[xPos]);
            }
            if (yMapping && residueMapping[yPos]) {
                residuesToHighlight.push(...residueMapping[yPos]);
            }

            // Reset all styles
            viewer.setStyle({}, { cartoon: { colorscheme: 'chain' } });

            // Highlight new residues
            if (residuesToHighlight.length > 0) {
                residuesToHighlight.forEach(res => {
                    viewer.setStyle(
                        { chain: res.chain, resi: res.resi },
                        { cartoon: { color: 'red' } }
                    );
                });
                viewer.render();
            }
        }
    }

    // Function to clear highlights
    function clearHighlights() {
        d3.selectAll('.residue')
            .classed('hover-highlight', false);

        // Reset protein structure highlights
        if (viewer) {
            viewer.setStyle({}, { cartoon: { colorscheme: 'chain' } });
            viewer.render();
        }
    }

    // Function to map matrix index to sequence position
    function mapIndexToSequence(index, sequencePositions) {
        for (let seqInfo of sequencePositions) {
            if (index >= seqInfo.start && index <= seqInfo.end) {
                return {
                    seqIndex: seqInfo.index,
                    resIndex: index - seqInfo.start
                };
            }
        }
        return null;
    }

    // Function to fetch PDB file and render the protein structure
    function fetchPDBAndRender(pdb_id, sequences) {
        const pdbUrl = `https://files.rcsb.org/download/${pdb_id}.pdb`;

        fetch(pdbUrl)
            .then(response => {
                if (!response.ok) {
                    throw new Error(`Network response was not ok for PDB ID ${pdb_id}`);
                }
                return response.text();
            })
            .then(pdbData => {
                parsePDBAndMapResidues(pdbData, sequences);
                renderProteinStructure(pdbData, pdb_id);
            })
            .catch(error => {
                console.error(`Error fetching PDB file for ${pdb_id}:`, error);
                d3.select('#protein-visualizer-container')
                    .append('div')
                    .attr('class', 'error')
                    .text(`Error loading PDB structure for ${pdb_id}: ${error.message}`);
            });
    }

    // Function to parse PDB and map residues
    function parsePDBAndMapResidues(pdbData, sequences) {
        residueMapping = {}; // Reset residue mapping

        const pdbLines = pdbData.split('\n');
        const chainSequences = {};
        const chainResidues = {};

        // Extract sequences from SEQRES records
        pdbLines.forEach(line => {
            if (line.startsWith('SEQRES')) {
                const chainID = line.substring(11, 12).trim();
                const resNames = line.substring(19).trim().split(/\s+/);
                if (!chainSequences[chainID]) {
                    chainSequences[chainID] = '';
                    chainResidues[chainID] = [];
                }
                resNames.forEach(resName => {
                    const resCode = threeLetterToOneLetter(resName);
                    chainSequences[chainID] += resCode;
                    chainResidues[chainID].push({ resi: null, resn: resName });
                });
            }
        });

        // Map residue numbers from ATOM records
        pdbLines.forEach(line => {
            if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
                const chainID = line.substring(21, 22).trim();
                const resSeq = parseInt(line.substring(22, 26).trim());
                const resName = line.substring(17, 20).trim();

                if (chainResidues[chainID]) {
                    // Find the first residue with matching resn and resi is null
                    const residue = chainResidues[chainID].find(res => res.resn === resName && res.resi === null);
                    if (residue) {
                        residue.resi = resSeq;
                    }
                }
            }
        });

        // Print sequences of each chain to console for debugging
        console.log("PDB Chain Sequences:");
        Object.keys(chainSequences).forEach(chainID => {
            const sequence = chainSequences[chainID];
            console.log(`Chain ${chainID}: ${sequence}`);
        });

        // Print input sequences to console for debugging
        console.log("Input Sequences:");
        sequences.forEach((seq, seqIndex) => {
            console.log(`Sequence ${seqIndex}: ${seq}`);
        });

        // Map input sequences to PDB chain sequences
        sequences.forEach((inputSeq, seqIndex) => {
            Object.keys(chainSequences).forEach(chainID => {
                const chainSeq = chainSequences[chainID];
                const alignmentIndex = simpleSequenceAlignment(chainSeq, inputSeq);

                if (alignmentIndex !== -1) {
                    console.log(`Sequence ${seqIndex} aligned with chain ${chainID} at position ${alignmentIndex}`);
                    // Map each residue position
                    for (let i = 0; i < inputSeq.length; i++) {
                        const globalIndex = sequences.slice(0, seqIndex).reduce((acc, seq) => acc + seq.length + '<+>'.length, 0) + i;
                        const chainResidue = chainResidues[chainID][alignmentIndex + i];

                        if (!residueMapping[globalIndex]) {
                            residueMapping[globalIndex] = [];
                        }

                        residueMapping[globalIndex].push({
                            chain: chainID,
                            resi: chainResidue.resi || (alignmentIndex + i + 1) // Assign resi if missing
                        });
                    }
                } else {
                    console.log(`Sequence ${seqIndex} did not align with chain ${chainID}`);
                }
            });
        });

        // For debugging: Print residueMapping
        console.log("Residue Mapping:", residueMapping);
    }

    // Simple sequence alignment function
    function simpleSequenceAlignment(chainSeq, inputSeq) {
        const maxMismatch = 5; // Allow up to 5 mismatches
        for (let i = 0; i <= chainSeq.length - inputSeq.length; i++) {
            let mismatches = 0;
            for (let j = 0; j < inputSeq.length; j++) {
                if (chainSeq[i + j] !== inputSeq[j]) {
                    mismatches++;
                    if (mismatches > maxMismatch) break;
                }
            }
            if (mismatches <= maxMismatch) {
                return i;
            }
        }
        return -1;
    }

    // Function to render protein structure using 3Dmol.js
    function renderProteinStructure(pdbData, pdb_id) {
        const element = document.getElementById('protein-visualizer');
        element.innerHTML = ''; // Clear previous visualization

        viewer = $3Dmol.createViewer(element, { defaultcolors: $3Dmol.rasmolElementColors });

        viewer.addModel(pdbData, 'pdb');

        // Collect unique chain IDs manually
        const chainIDs = new Set();
        viewer.getModel().atoms.forEach(atom => {
            if (atom.chain && !chainIDs.has(atom.chain)) {
                chainIDs.add(atom.chain);
            }
        });

        // Print chain IDs and their sequences for debugging
        console.log("Unique Chain IDs:", Array.from(chainIDs));

        // Assign unique colors to each chain
        const colorScale = d3.scaleOrdinal(d3.schemeCategory10);
        let chainIndex = 0;
        chainIDs.forEach(chainID => {
            viewer.setStyle(
                { chain: chainID },
                { cartoon: { color: colorScale(chainIndex) } }
            );
            chainIndex++;
        });

        // Highlight residues based on residueMapping
        // (Residues will be highlighted on hover events)

        viewer.zoomTo();
        viewer.render();
    }

    // Helper function to convert three-letter amino acid codes to one-letter codes
    function threeLetterToOneLetter(threeLetter) {
        const aaTable = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
            'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
            'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
            'XLE': 'J', 'XAA': 'X', 'UNK': 'X', 'MSE': 'M'
        };
        return aaTable[threeLetter.toUpperCase()] || 'X';
    }

});
