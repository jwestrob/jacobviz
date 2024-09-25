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
    let sequenceOffsets = {};

    // Mapping of chainID to color
    let chainColorMap = {};

    // Map of 3-letter codes to 1-letter codes
    const threeToOne = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    };

    // Variables for performance optimization
    let isZooming = false;
    let mouseoverTimeout = null;
    let renderTimeout = null;
    let previouslyHighlightedResidues = [];

    // Variables for CIF data and atoms
    let cifData = null;
    let cifAtoms = null;

    // Define a list of standard amino acid residue names to exclude waters and other unwanted residues
    const standardResidues = new Set([
        'ALA', 'CYS', 'ASP', 'GLU', 'PHE',
        'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
        'MET', 'ASN', 'PRO', 'GLN', 'ARG',
        'SER', 'THR', 'VAL', 'TRP', 'TYR',
        // Add nucleic acids if necessary
        'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'T', 'U'
    ]);

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
        const middleContainer = document.getElementById('middle-container');
        const containerWidth = middleContainer.clientWidth;
        const containerHeight = middleContainer.clientHeight;

        // Since CSS handles the layout, we can set heatmapWidth and heatmapHeight based on containerWidth
        // Assuming heatmap-container and protein-visualizer-container each take ~48% of the width
        heatmapWidth = middleContainer.clientWidth * 0.48 - 20; // 20 accounts for margin
        heatmapHeight = middleContainer.clientHeight - 20; // Adjust as needed

        // Log dimensions for debugging
        console.log('Middle Container Width:', containerWidth);
        console.log('Heatmap Width:', heatmapWidth);
        console.log('Heatmap Height:', heatmapHeight);

        // Create scales
        const xScale = d3.scaleBand()
            .domain(d3.range(matrixSize))
            .range([0, heatmapWidth])
            .padding(0.01);

        const yScale = d3.scaleBand()
            .domain(d3.range(matrixSize))
            .range([0, heatmapHeight]) // Adjusted Y-axis range
            .padding(0.01);

        const colorScale = d3.scaleSequential(d3.interpolateViridis)
            .domain([d3.min(matrixData.flat()), d3.max(matrixData.flat())]);

        // Render sequences
        renderSequences(sequences, sequencePositions);

        // Render Jacobian heatmap
        renderHeatmap(matrixData, xScale, yScale, colorScale, sequencePositions);

        // Render static Jacobian contact map at the bottom right
        renderStaticJacobianContactMap(matrixData);

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

    // Function to render the static Jacobian contact map
    function renderStaticJacobianContactMap(matrixData) {
        const container = d3.select('#jacobian-contact-map-container');
        container.html(''); // Clear existing content

        // Add a title
        container.append('h2').text('Jacobian Contact Map');

        // Dimensions
        const margin = { top: 50, right: 50, bottom: 50, left: 50 };
        const width = heatmapWidth + margin.left + margin.right;
        const height = heatmapHeight + margin.top + margin.bottom;
        const matrixSize = matrixData.length;

        // Create scales
        const xScale = d3.scaleBand()
            .domain(d3.range(matrixSize))
            .range([0, heatmapWidth]);

        const yScale = d3.scaleBand()
            .domain(d3.range(matrixSize))
            .range([0, heatmapHeight]);

        // Create SVG element
        const svg = container.append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // Flatten matrix data for easier processing
        const flatData = [];
        for (let y = 0; y < matrixSize; y++) {
            for (let x = 0; x < matrixSize; x++) {
                flatData.push({
                    value: matrixData[y][x],
                    x: x,
                    y: y
                });
            }
        }

        // Define color scale
        const colorScaleStatic = d3.scaleSequential(d3.interpolateViridis)
            .domain([d3.min(matrixData.flat()), d3.max(matrixData.flat())]);

        // Render heatmap cells
        svg.selectAll('rect')
            .data(flatData)
            .enter()
            .append('rect')
            .attr('x', d => xScale(d.x))
            .attr('y', d => yScale(d.y))
            .attr('width', xScale.bandwidth())
            .attr('height', yScale.bandwidth())
            .style('fill', d => colorScaleStatic(d.value));

        // Add x-axis
        const xAxis = d3.axisBottom(xScale)
            .tickValues(xScale.domain().filter((d, i) => !(i % 10))) // Every 10 units
            .tickFormat(d => d);

        svg.append('g')
            .attr('class', 'x axis')
            .attr('transform', `translate(0,${heatmapHeight})`)
            .call(xAxis)
            .selectAll("text")
            .attr("transform", "rotate(-90)")
            .attr("dx", "-0.8em")
            .attr("dy", "-0.5em")
            .style("text-anchor", "end");

        // Add y-axis
        const yAxis = d3.axisLeft(yScale)
            .tickValues(yScale.domain().filter((d, i) => !(i % 10))) // Every 10 units
            .tickFormat(d => d);

        svg.append('g')
            .attr('class', 'y axis')
            .call(yAxis);

        // Add axis labels (only x-axis label as per your request)
        svg.append('text')
            .attr('x', heatmapWidth / 2)
            .attr('y', heatmapHeight + margin.bottom - 10)
            .attr('text-anchor', 'middle')
            .text('Residue Index');
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

        sequences.forEach((seq, seqIndex) => {
            const seqDiv = sequenceContainer.append('div')
                .attr('id', `sequence-${seqIndex}`)
                .attr('class', 'sequence');

            // Add sequence label
            seqDiv.append('div')
                .attr('class', 'sequence-label')
                .text(`Sequence ${seqIndex + 1}:`)
                .style('margin-bottom', '5px');

            // Create a div for the sequence content
            const seqContentDiv = seqDiv.append('div')
                .attr('class', 'sequence-content');

            // Split the sequence into lines of 150 characters
            const lineLength = 150;
            for (let i = 0; i < seq.length; i += lineLength) {
                const lineSeq = seq.substring(i, i + lineLength);

                // Create a span for the line
                const lineSpan = seqContentDiv.append('div')
                    .attr('class', 'sequence-line');

                lineSeq.split('').forEach((residue, resIndex) => {
                    const globalResIndex = i + resIndex;
                    lineSpan.append('span')
                        .attr('id', `seq${seqIndex}-res${globalResIndex}`)
                        .attr('class', 'residue')
                        .text(residue)
                        .on('click', function () {
                            const element = d3.select(this);
                            element.classed('manual-highlight', !element.classed('manual-highlight'));
                        })
                        .on('mouseover', function () {
                            const index = sequencePositions[seqIndex].start + globalResIndex;
                            highlightResidues(index, index, sequencePositions);
                        })
                        .on('mouseout', function () {
                            clearHighlights();
                        });
                });
            }
        });
    }

    // Function to render heatmap using SVG
    function renderHeatmap(matrixData, xScale, yScale, colorScale, sequencePositions) {
        const heatmapContainer = d3.select('#heatmap-container');
        heatmapContainer.html(''); // Clear existing heatmap

        const margin = { top: 50, right: 50, bottom: 50, left: 50 };

        const svg = heatmapContainer.append('svg')
            .attr('width', heatmapWidth + margin.left + margin.right)
            .attr('height', heatmapHeight + margin.top + margin.bottom)
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
        const cells = heatmapGroup.selectAll('rect')
            .data(flatData)
            .enter()
            .append('rect')
            .attr('x', d => xScale(d.x))
            .attr('y', d => yScale(d.y))
            .attr('width', xScale.bandwidth())
            .attr('height', yScale.bandwidth())
            .style('fill', d => colorScale(d.value))
            .on('mouseover', function (event, d) {
                if (isZooming) return; // Skip if zooming
                clearTimeout(mouseoverTimeout);
                mouseoverTimeout = setTimeout(() => {
                    handleMouseOver(event, d, sequencePositions);
                }, 100); // Debounce delay (in milliseconds)
            })
            .on('mouseout', function () {
                clearTimeout(mouseoverTimeout);
                handleMouseOut();
            });

        // Remove axes to prevent overlap when zoomed

        // Zoom function without axis updates
        const zoom = d3.zoom()
            .scaleExtent([1, 10])
            .on('start', zoomStarted)
            .on('zoom', zoomed)
            .on('end', zoomEnded);

        svg.call(zoom);

        function zoomStarted() {
            isZooming = true;
        }

        function zoomed(event) {
            const transform = event.transform;
            heatmapGroup.attr('transform', transform);
        }

        function zoomEnded() {
            isZooming = false;
        }
    }

    // Function to handle mouse over event
    function handleMouseOver(event, d, sequencePositions) {
        // Highlight corresponding residues
        highlightResidues(d.x, d.y, sequencePositions);

        // Get adjusted indices
        const xMapping = mapIndexToSequence(d.x, sequencePositions);
        const yMapping = mapIndexToSequence(d.y, sequencePositions);

        const adjustedXIndex = xMapping ? xMapping.adjustedIndex + 1 : d.x + 1;
        const adjustedYIndex = yMapping ? yMapping.adjustedIndex + 1 : d.y + 1;

        // Display tooltip at the mouse position
        const tooltip = d3.select('#tooltip');
        tooltip.classed('hidden', false)
            .classed('visible', true)
            .html(`Value: ${d.value.toFixed(4)}<br>X: ${adjustedXIndex}<br>Y: ${adjustedYIndex}`)
            .style('left', (event.pageX + 15) + 'px')
            .style('top', (event.pageY - 28) + 'px');
    }

    // Function to handle mouse out event
    function handleMouseOut() {
        clearHighlights();
        d3.select('#tooltip')
            .classed('hidden', true)
            .classed('visible', false);

        // Clear residue info
        d3.select('#residue-info').html('');
    }

    // Function to highlight residues
    function highlightResidues(xPos, yPos, sequencePositions) {
        const newHighlightedResidues = [];

        const xMapping = mapIndexToSequence(xPos, sequencePositions);
        const yMapping = mapIndexToSequence(yPos, sequencePositions);

        let residueInfoText = '';

        if (xMapping) {
            const residueLetter = sequencePositions[xMapping.seqIndex].sequence[xMapping.resIndex];
            d3.select(`#seq${xMapping.seqIndex}-res${xMapping.resIndex}`)
                .classed('hover-highlight', true);
            residueInfoText += `Residue X: Index ${xMapping.adjustedIndex + 1}, Letter ${residueLetter}<br>`;
        }
        if (yMapping) {
            const residueLetter = sequencePositions[yMapping.seqIndex].sequence[yMapping.resIndex];
            d3.select(`#seq${yMapping.seqIndex}-res${yMapping.resIndex}`)
                .classed('hover-highlight', true);
            residueInfoText += `Residue Y: Index ${yMapping.adjustedIndex + 1}, Letter ${residueLetter}<br>`;
        }

        // Update residue info display
        d3.select('#residue-info').html(residueInfoText);

        // Highlight residues in the protein structure
        if (viewer && residueMapping) {
            if (xMapping && residueMapping[xPos]) {
                newHighlightedResidues.push(residueMapping[xPos]);
            }
            if (yMapping && residueMapping[yPos]) {
                newHighlightedResidues.push(residueMapping[yPos]);
            }

            // Remove highlighting from previously highlighted residues
            previouslyHighlightedResidues.forEach(res => {
                // Reset to the chain's original color
                const chainColor = chainColorMap[res.chain];
                if (chainColor) {
                    viewer.setStyle(
                        { chain: res.chain, resi: res.resi },
                        { cartoon: { color: chainColor } }
                    );
                } else {
                    // If no chain color found, reset to default
                    viewer.setStyle(
                        { chain: res.chain, resi: res.resi },
                        { cartoon: { color: 'grey' } }
                    );
                }
            });

            // Apply highlighting to new residues
            newHighlightedResidues.forEach(res => {
                viewer.setStyle(
                    { chain: res.chain, resi: res.resi },
                    { cartoon: { color: 'red' } }
                );
            });

            // Update the list of currently highlighted residues
            previouslyHighlightedResidues = newHighlightedResidues;

            // Render only if necessary
            viewer.render();
        }
    }

    // Function to clear highlights
    function clearHighlights() {
        // Clear highlights in the sequence
        d3.selectAll('.residue.hover-highlight').classed('hover-highlight', false);

        // Reset styles in the 3D structure
        if (viewer) {
            viewer.setStyle({}, {}); // Clear all styles
            // Reapply default styles to all chains
            for (const chainID in chainColorMap) {
                viewer.setStyle(
                    { chain: chainID },
                    { cartoon: { color: chainColorMap[chainID] } }
                );
            }
            viewer.render();
        }
    }

    // Function to map matrix index to sequence position
    function mapIndexToSequence(index, sequencePositions) {
        for (let seqInfo of sequencePositions) {
            if (index >= seqInfo.start && index <= seqInfo.end) {
                const seqIndex = seqInfo.index;
                const resIndex = index - seqInfo.start;
                const offset = sequenceOffsets[seqIndex] || 0;
                const adjustedIndex = resIndex - offset;
                return {
                    seqIndex: seqIndex,
                    resIndex: resIndex,
                    adjustedIndex: adjustedIndex
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
                // Offload parsing to a Web Worker
                const worker = new Worker('pdbWorker.js');
                worker.postMessage({ pdbData, sequences });

                worker.onmessage = function (e) {
                    residueMapping = e.data.residueMapping;
                    sequenceOffsets = e.data.sequenceOffsets; // Store sequence offsets
                    renderProteinStructure(pdbData, pdb_id);
                    worker.terminate(); // Terminate the worker after completion

                    // After rendering the protein structure, fetch and prepare the CIF contact map controls
                    fetchCIFAndPrepareContactMap(pdb_id);
                };

                worker.onerror = function (error) {
                    console.error(`Worker error:`, error);
                    worker.terminate();
                };
            })
            .catch(error => {
                console.error(`Error fetching PDB file for ${pdb_id}:`, error);
                d3.select('#protein-visualizer-container')
                    .append('div')
                    .attr('class', 'error')
                    .text(`Error loading PDB structure for ${pdb_id}: ${error.message}`);
            });
    }

    // Function to render protein structure using 3Dmol.js
    function renderProteinStructure(pdbData, pdb_id) {
        const element = document.getElementById('protein-visualizer');
        element.innerHTML = ''; // Clear previous visualization

        viewer = $3Dmol.createViewer(element, {
            defaultcolors: $3Dmol.rasmolElementColors,
            width: element.clientWidth,
            height: element.clientHeight
        });

        // Add model with doAssembly option set to true
        viewer.addModel(pdbData, 'pdb', { doAssembly: true });

        // Collect unique chain IDs manually
        const chainIDs = new Set();
        viewer.getModel().atoms.forEach(atom => {
            if (atom.chain && !chainIDs.has(atom.chain)) {
                chainIDs.add(atom.chain);
            }
        });

        // Assign unique colors to each chain and store in chainColorMap
        chainColorMap = {}; // Reset chainColorMap
        const colorScale = d3.scaleOrdinal(d3.schemeCategory10);
        let chainIndex = 0;
        chainIDs.forEach(chainID => {
            const color = colorScale(chainIndex);
            chainColorMap[chainID] = color;
            viewer.setStyle(
                { chain: chainID, resn: Object.keys(threeToOne) }, // Only standard amino acids
                { cartoon: { color: color } }
            );
            chainIndex++;
        });

        viewer.zoomTo();
        viewer.render();

        // Handle window resize to adjust viewer size
        window.addEventListener('resize', function () {
            if (viewer) {
                viewer.resize();
                viewer.render();
            }
        });
    }
    
    function getResidueInfo(i, j, sequencePositions) {
        const seqInfoI = mapIndexToSequence(i, sequencePositions);
        const seqInfoJ = mapIndexToSequence(j, sequencePositions);

        const residueI = residueMapping[i];
        const residueJ = residueMapping[j];

        const adjustedIndexI = seqInfoI ? seqInfoI.adjustedIndex + 1 : i + 1;
        const adjustedIndexJ = seqInfoJ ? seqInfoJ.adjustedIndex + 1 : j + 1;

        const resInfoI = residueI ? `Chain ${residueI.chain}, Residue ${residueI.resi}` : 'N/A';
        const resInfoJ = residueJ ? `Chain ${residueJ.chain}, Residue ${residueJ.resi}` : 'N/A';

        return `Position ${adjustedIndexI} (${resInfoI}) ↔ Position ${adjustedIndexJ} (${resInfoJ})`;
    }

    // Function to fetch CIF file and prepare for contact map generation
    function fetchCIFAndPrepareContactMap(pdb_id) {
        const cifUrl = `https://files.rcsb.org/download/${pdb_id}.cif`;

        fetch(cifUrl)
            .then(response => {
                if (!response.ok) {
                    throw new Error(`Network response was not ok for CIF file of PDB ID ${pdb_id}`);
                }
                return response.text();
            })
            .then(cifText => {
                cifData = cifText;
                prepareContactMapControls(cifData);
            })
            .catch(error => {
                console.error(`Error fetching CIF file for ${pdb_id}:`, error);
                d3.select('#contact-map-container')
                    .append('div')
                    .attr('class', 'error')
                    .text(`Error loading CIF file for ${pdb_id}: ${error.message}`);
            });
    }

    // Function to prepare contact map controls and populate chain options
    function prepareContactMapControls(cifData) {
        // Create a new GLModel instance with doAssembly option
        const model = new $3Dmol.GLModel(0, { doAssembly: true });

        // Add the CIF data to the model
        model.addMolData(cifData, 'cif');

        // Get the atoms
        cifAtoms = model.atoms;

        // Extract unique chain IDs
        const chainIDs = new Set();
        cifAtoms.forEach(atom => {
            if (atom.chain) {
                chainIDs.add(atom.chain);
            }
        });

        // Populate the chain-select dropdown
        const chainSelect = document.getElementById('chain-select');
        chainSelect.innerHTML = ''; // Clear previous options

        chainIDs.forEach(chainID => {
            const option = document.createElement('option');
            option.value = chainID;
            option.text = chainID;
            chainSelect.appendChild(option);
        });

        // By default, select all chains
        for (let i = 0; i < chainSelect.options.length; i++) {
            chainSelect.options[i].selected = true;
        }

        // Enable the "Generate Contact Map" button
        const generateButton = document.getElementById('generate-contact-map-button');
        generateButton.disabled = false;

        // Remove previous event listeners to prevent multiple bindings
        generateButton.replaceWith(generateButton.cloneNode(true));
        const newGenerateButton = document.getElementById('generate-contact-map-button');

        // Set up event listener for the button
        newGenerateButton.addEventListener('click', function () {
            const selectedChains = Array.from(chainSelect.selectedOptions).map(option => option.value);
            const distanceThresholdInput = document.getElementById('distance-threshold');
            const distanceThreshold = parseFloat(distanceThresholdInput.value);

            // Debug: Log selected chains and distance threshold
            console.log(`Selected Chains: ${selectedChains}`);
            console.log(`Distance Threshold: ${distanceThreshold} Å`);

            generateCIFContactMap(cifData, selectedChains, distanceThreshold);
        });
    }

    // Function to generate contact map from CIF data
    function generateCIFContactMap(cifData, selectedChains, distanceThreshold) {
        // Use the stored cifAtoms
        const atoms = cifAtoms;

        // Filter atoms by selected chains and exclude waters ('HOH')
        const filteredAtoms = atoms.filter(atom => {
            return selectedChains.includes(atom.chain) && standardResidues.has(atom.resn);
        });

        // Debug: Log number of atoms after filtering
        console.log(`Number of atoms after filtering by chains (${selectedChains.join(', ')}), excluding waters: ${filteredAtoms.length}`);

        // Map to store residues
        const residuesMap = new Map();

        filteredAtoms.forEach(atom => {
            if (atom.elem === 'H') return; // Skip hydrogen atoms

            const chain = atom.chain || '';
            const resi = atom.resi;
            const resn = atom.resn;

            // Exclude residues not in standardResidues
            if (!standardResidues.has(resn)) {
                return;
            }

            const residueKey = `${chain}_${resi}`;

            if (!residuesMap.has(residueKey)) {
                residuesMap.set(residueKey, {
                    chain: chain,
                    resi: resi,
                    resn: resn,
                    atoms: []
                });
            }

            residuesMap.get(residueKey).atoms.push(atom);
        });

        // Debug: Log number of residues after grouping
        console.log(`Number of residues after grouping (excluding waters): ${residuesMap.size}`);

        // Now, for each residue, calculate the positions of all heavy atoms
        const residues = [];
        const residueIndices = []; // Map to original indices
        let index = 0;

        residuesMap.forEach(residue => {
            const heavyAtoms = residue.atoms.filter(atom => atom.elem !== 'H');
            if (heavyAtoms.length > 0) {
                residues.push({
                    chain: residue.chain,
                    resi: residue.resi,
                    resn: residue.resn,
                    atoms: heavyAtoms
                });
                residueIndices.push(index);
                index++;
            }
        });

        // Debug: Log residues details
        console.log(`Residues details (excluding waters):`);
        residues.forEach((residue, resIndex) => {
            console.log(`Residue ${resIndex}: Chain ${residue.chain}, Residue ${residue.resi}, Residue Name ${residue.resn}`);
        });

        // Extract coordinates of all heavy atoms
        const coords = [];
        const residueAtomIndices = [];
        residues.forEach((residue, resIndex) => {
            residue.atoms.forEach(atom => {
                coords.push([atom.x, atom.y, atom.z]);
                residueAtomIndices.push(resIndex);
            });
        });

        // Initialize contact map with NaN
        const totalResidues = residues.length;
        const contactMap = new Array(totalResidues);
        for (let i = 0; i < totalResidues; i++) {
            contactMap[i] = new Array(totalResidues).fill(NaN);
        }

        if (coords.length > 0) {
            const coordsArray = coords; // Array of [x, y, z]
            const numAtoms = coordsArray.length;

            // Calculate pairwise distances between all atoms
            for (let i = 0; i < numAtoms; i++) {
                const xi = coordsArray[i][0];
                const yi = coordsArray[i][1];
                const zi = coordsArray[i][2];
                const resIndexI = residueAtomIndices[i];

                for (let j = i; j < numAtoms; j++) {
                    const xj = coordsArray[j][0];
                    const yj = coordsArray[j][1];
                    const zj = coordsArray[j][2];
                    const resIndexJ = residueAtomIndices[j];

                    const dx = xi - xj;
                    const dy = yi - yj;
                    const dz = zi - zj;
                    const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

                    if (distance <= distanceThreshold) {
                        contactMap[resIndexI][resIndexJ] = 1;
                        contactMap[resIndexJ][resIndexI] = 1;
                    }
                }
            }

            // Set observed residues to 0 where there is no contact
            for (let i = 0; i < totalResidues; i++) {
                for (let j = 0; j < totalResidues; j++) {
                    if (isNaN(contactMap[i][j])) {
                        contactMap[i][j] = 0;
                    }
                }
            }
        }

        // Debug: Log the size of the contact map
        console.log(`Generated contact map size: ${contactMap.length} x ${contactMap.length}`);

        // Debug: Compare with expected residue count
        // Replace with dynamic value if available
        const expectedResidues = 359; // For PDB ID 3kzk
        console.log(`Expected number of residues: ${expectedResidues}`);
        console.log(`Actual number of residues in contact map: ${totalResidues}`);

        // Log the sequences being used
        const sequencesUsed = residues.map(res => res.resn);
        console.log(`Sequences used for contact map: ${sequencesUsed.join('')}`);

        // Now we can render the contact map
        renderContactMap(contactMap, residues);
    }

    // Function to render the contact map
    function renderContactMap(contactMap, residues) {
        const container = d3.select('#cif-contact-map-container');
        container.html(''); // Clear existing content

        // Add a title
        container.append('h2').text('Contact Map from CIF File');

        // Dimensions
        const margin = { top: 50, right: 50, bottom: 50, left: 50 };
        const width = heatmapWidth + margin.left + margin.right;
        const height = heatmapHeight + margin.top + margin.bottom;
        const totalResidues = contactMap.length;

        // Create scales
        const xScale = d3.scaleBand()
            .domain(d3.range(totalResidues))
            .range([0, heatmapWidth]);

        const yScale = d3.scaleBand()
            .domain(d3.range(totalResidues))
            .range([0, heatmapHeight]);

        // Create SVG element
        const svg = container.append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // Flatten contact map data
        const flatData = [];
        for (let i = 0; i < totalResidues; i++) {
            for (let j = 0; j < totalResidues; j++) {
                flatData.push({
                    x: i,
                    y: j,
                    value: contactMap[i][j]
                });
            }
        }

        // Define color scale
        const colorScale = d3.scaleSequential(d3.interpolateBlues)
            .domain([0, 1]);

        // Render heatmap cells
        svg.selectAll('rect')
            .data(flatData)
            .enter()
            .append('rect')
            .attr('x', d => xScale(d.x))
            .attr('y', d => yScale(d.y))
            .attr('width', xScale.bandwidth())
            .attr('height', yScale.bandwidth())
            .style('fill', d => colorScale(d.value));

        // Add axes with residue indices
        const xAxis = d3.axisBottom(xScale)
            .tickValues(xScale.domain().filter((d, i) => !(i % 10))) // Every 10 residues
            .tickFormat(d => residues[d].resi);

        const yAxis = d3.axisLeft(yScale)
            .tickValues(yScale.domain().filter((d, i) => !(i % 10))) // Every 10 residues
            .tickFormat(d => residues[d].resi);

        svg.append('g')
            .attr('class', 'x axis')
            .attr('transform', `translate(0,${heatmapHeight})`)
            .call(xAxis)
            .selectAll("text")
            .attr("transform", "rotate(-90)")
            .attr("dx", "-0.8em")
            .attr("dy", "-0.5em")
            .style("text-anchor", "end");

        svg.append('g')
            .attr('class', 'y axis')
            .call(yAxis);

        // Add axis labels (removed y-axis label as per your request)
        svg.append('text')
            .attr('x', heatmapWidth / 2)
            .attr('y', heatmapHeight + margin.bottom - 10)
            .attr('text-anchor', 'middle')
            .text('Residue Index');
    }


});
