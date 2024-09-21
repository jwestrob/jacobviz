// script.js

// Ensure the script runs after the page has loaded
document.addEventListener('DOMContentLoaded', function () {
    // API endpoint to get the list of JSON data files
    const dataFilesEndpoint = '/api/data-files';

    // Heatmap dimensions
    const heatmapWidth = 700;
    const heatmapHeight = 700;

    // Initialize a cache object to store fetched JSON data
    const dataCache = {};

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
        renderSequences(sequences);

        // Render heatmap
        renderHeatmap(matrixData, xScale, yScale, colorScale, sequencePositions);
    }

    // Function to parse sequences
    function parseSequences(sequenceInput) {
        return sequenceInput.split('<+>').filter(seq => seq.length > 0);
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
    function renderSequences(sequences) {
        const sequenceContainer = d3.select('#sequences-container');

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
                    });
            });
        });
    }

    // Function to render heatmap
    function renderHeatmap(matrixData, xScale, yScale, colorScale, sequencePositions) {
        const heatmapWidth = 700;
        const heatmapHeight = 700;
        const margin = { top: 20, right: 20, bottom: 20, left: 20 };

        const svg = d3.select('#heatmap-container')
            .append('svg')
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
        heatmapGroup.selectAll()
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
    }

    // Function to handle mouse over event
    function handleMouseOver(event, d, sequencePositions, xScale, yScale, svg) {
        // Get the current transform
        const transform = d3.zoomTransform(svg.node());

        // Adjust the position based on the transform
        const xPosition = transform.applyX(xScale(d.x) + xScale.bandwidth() / 2);
        const yPosition = transform.applyY(yScale(d.y) + yScale.bandwidth() / 2);

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
    }

    // Function to clear highlights
    function clearHighlights() {
        d3.selectAll('.residue')
            .classed('hover-highlight', false);
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

});
