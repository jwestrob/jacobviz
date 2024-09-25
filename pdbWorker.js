// pdbWorker.js

self.onmessage = function (e) {
    const { pdbData, sequences } = e.data;
    const residueMapping = parsePDBAndMapResidues(pdbData, sequences);
    self.postMessage({ residueMapping });
};

function parsePDBAndMapResidues(pdbData, sequences) {
    let residueMapping = {}; // Initialize residue mapping

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

    return residueMapping;
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
