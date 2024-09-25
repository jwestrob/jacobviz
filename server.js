const express = require('express');
const fs = require('fs');
const path = require('path');
const cors = require('cors');

const app = express();
const PORT = 3000;

// Enable CORS if accessing from a different domain
app.use(cors());

// Serve static files from the root directory
app.use(express.static(path.join(__dirname)));

// API endpoint to list all JSON files in the data/ directory
app.get('/api/data-files', (req, res) => {
    const dataDir = path.join(__dirname, 'data');
    fs.readdir(dataDir, (err, files) => {
        if (err) {
            console.error('Error reading data directory:', err);
            return res.status(500).json({ error: 'Unable to list data files.' });
        }
        // Filter out JSON files
        const jsonFiles = files.filter(file => path.extname(file).toLowerCase() === '.json');
        res.json({ data_files: jsonFiles });
    });
});

// Start the server
app.listen(PORT, () => {
    console.log(`Server is running at http://localhost:${PORT}`);
});
