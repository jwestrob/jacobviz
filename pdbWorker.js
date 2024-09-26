// pdbWorker.js

// List of standard amino acid residue names
const standardAminoAcids = new Set([
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
]);

// Map of 3-letter codes to 1-letter codes
const threeToOne = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
};

onmessage = function (e) {
    const pdbData = e.data.pdbData;
    const sequences = e.data.sequences;
    const residueMapping = {};
    const sequenceOffsets = {}; // To store offsets for each sequence

    // Parse the PDB file and extract sequences for each chain
    const chainResidues = {};
    const lines = pdbData.split('\n');
    lines.forEach(line => {
        if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
            const resName = line.substring(17, 20).trim();
            const chainID = line.substring(21, 22);
            const resSeq = parseInt(line.substring(22, 26).trim());

            // Skip heteroatoms and non-standard amino acids
            if (!standardAminoAcids.has(resName)) return;

            // Initialize chain if not present
            if (!chainResidues[chainID]) {
                chainResidues[chainID] = [];
            }

            // Avoid duplicates
            const lastResidue = chainResidues[chainID][chainResidues[chainID].length - 1];
            if (lastResidue && lastResidue.resSeq === resSeq) return;

            chainResidues[chainID].push({
                resName: resName,
                resSeq: resSeq
            });
        }
    });

    // Build sequences for each chain
    const chainSequences = {};
    for (const chainID in chainResidues) {
        const residues = chainResidues[chainID];
        let sequence = '';
        residues.forEach(res => {
            sequence += threeToOne[res.resName] || 'X';
        });
        chainSequences[chainID] = sequence;
    }

    // Implement Needleman-Wunsch algorithm for sequence alignment
    function needlemanWunsch(seq1, seq2) {
        const matchScore = 1;
        const mismatchScore = -1;
        const gapPenalty = -1;

        const n = seq1.length;
        const m = seq2.length;

        const scoreMatrix = [];
        const traceMatrix = [];

        // Initialize matrices
        for (let i = 0; i <= n; i++) {
            scoreMatrix[i] = [];
            traceMatrix[i] = [];
            for (let j = 0; j <= m; j++) {
                if (i === 0) {
                    scoreMatrix[i][j] = j * gapPenalty;
                    traceMatrix[i][j] = 'left';
                } else if (j === 0) {
                    scoreMatrix[i][j] = i * gapPenalty;
                    traceMatrix[i][j] = 'up';
                } else {
                    scoreMatrix[i][j] = 0;
                    traceMatrix[i][j] = '';
                }
            }
        }

        // Fill matrices
        for (let i = 1; i <= n; i++) {
            for (let j = 1; j <= m; j++) {
                const match = scoreMatrix[i - 1][j - 1] + (seq1[i - 1] === seq2[j - 1] ? matchScore : mismatchScore);
                const deleteGap = scoreMatrix[i - 1][j] + gapPenalty;
                const insertGap = scoreMatrix[i][j - 1] + gapPenalty;

                const maxScore = Math.max(match, deleteGap, insertGap);
                scoreMatrix[i][j] = maxScore;

                if (maxScore === match) {
                    traceMatrix[i][j] = 'diag';
                } else if (maxScore === deleteGap) {
                    traceMatrix[i][j] = 'up';
                } else {
                    traceMatrix[i][j] = 'left';
                }
            }
        }

        // Traceback
        let align1 = '';
        let align2 = '';
        let i = n;
        let j = m;

        while (i > 0 || j > 0) {
            if (traceMatrix[i][j] === 'diag') {
                align1 = seq1[i - 1] + align1;
                align2 = seq2[j - 1] + align2;
                i--;
                j--;
            } else if (traceMatrix[i][j] === 'up') {
                align1 = seq1[i - 1] + align1;
                align2 = '-' + align2;
                i--;
            } else {
                align1 = '-' + align1;
                align2 = seq2[j - 1] + align2;
                j--;
            }
        }

        return { align1, align2 };
    }

    // Map sequence indices to PDB residues using alignment
    let globalSeqIndex = 0;
    sequences.forEach((seq, seqIdx) => {
        let matched = false;
        let bestScore = -Infinity;
        let bestAlignment = null;
        let bestChainID = null;
        let bestResidues = null;

        for (const chainID in chainSequences) {
            const chainSeq = chainSequences[chainID];
            const { align1, align2 } = needlemanWunsch(seq, chainSeq);

            // Calculate alignment score
            let score = 0;
            for (let k = 0; k < align1.length; k++) {
                if (align1[k] === align2[k]) {
                    score += 1;
                } else if (align1[k] === '-' || align2[k] === '-') {
                    score -= 1;
                } else {
                    score -= 1;
                }
            }

            if (score > bestScore) {
                bestScore = score;
                bestAlignment = { align1, align2 };
                bestChainID = chainID;
                bestResidues = chainResidues[chainID];
            }
        }

        if (bestAlignment) {
            const { align1, align2 } = bestAlignment;
            let seqPos = 0;
            let pdbPos = 0;
            let seqOffset = null;

            for (let k = 0; k < align1.length; k++) {
                const aa1 = align1[k];
                const aa2 = align2[k];

                if (aa1 !== '-' && aa2 !== '-') {
                    if (seqOffset === null) {
                        seqOffset = seqPos - pdbPos;
                        sequenceOffsets[seqIdx] = seqOffset;
                    }
                residueMapping[globalSeqIndex + seqPos] = {
                    chain: bestChainID,
                    resi: bestResidues[pdbPos].resSeq,
                    adjustedIndex: seqPos - seqOffset,
                    residueLetter: align1[k]
                };
                    seqPos++;
                    pdbPos++;
                } else if (aa1 !== '-' && aa2 === '-') {
                    seqPos++;
                } else if (aa1 === '-' && aa2 !== '-') {
                    pdbPos++;
                }
            }
            matched = true;
        }

        if (!matched) {
            console.warn(`Sequence not matched with PDB chains: ${seq}`);
        }

        globalSeqIndex += seq.length;
    });

    postMessage({ residueMapping: residueMapping, sequenceOffsets: sequenceOffsets });
};
