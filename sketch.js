// --- Core Microbial Strategy Simulation (v29 - Stripes Gradient & Simplified Gal Cost) ---
// Features:
// - Simplified Galactose Cost: h_max now represents max relative rate for p=0 cell.
// - Added Stripe Gradients (Horizontal/Vertical) to test sharp transitions.
// - Multiple static nutrient gradient options with fine-grained control.
// - Evolving cell strategy 'p' & p-dependent Galactose cost scaling h_max.
// - Colonization Model: Viability Threshold, Collective Pressure, Probabilistic Success, Parent Limitation.
// - User interaction: Click to plant, Spacebar Start/Pause, Gradient Overlay Toggle.
// Optimizations: Frontier Opt, Buffers, Contextual UI.
// --- NO CELL DEATH in this version ---
// ---

// --- Simulation Settings ---
let maxInitRadius = 50;
let hexSize = 5;
const epsilon = 1e-9;
let drawInterval = 2;

// --- MODEL PARAMETERS (Defaults, adjustable via UI) ---
let PARAMS = {
  // Cell & Colonization
  initialP: 0.5,
  mutationProbability: 0.01,
  K_G: 0.5, K_H: 0.5,
  h_max: 0.85, // *** Represents MAX Gal rate relative to Glc for a p=0 cell ***
  // REMOVED: galactoseEfficiencyFactor
  minColonizationFitness: 0.05,
  colonizationPressureFactor: 0.8,
  maxColonizationProbability: 0.95,

  // Gradient Core
  gradientMode: "stripes_v", // Default to stripes

  // Smooth Radial/Linear Controls
  gradientCenter: 0.25,
  gradientTransitionWidth: 0.5,

  // Gaussian Source Controls
  gaussianSourceQ: 0, gaussianSourceR: 0, gaussianWidth: 15,
  gaussianPeakG: 1.0, gaussianPeakH: 0.0, gaussianBaseG: 0.0, gaussianBaseH: 0.2,

  // Noisy Gradient Controls
  noiseScale: 0.1, noiseIntensityG: 0.3, noiseIntensityH: 0.3,
  noiseBaseGradient: "uniform_mix",

  // *** NEW: Stripe Controls ***
  stripeWidth: 40, // Width of stripes in pixels (adjust based on hexSize/visuals)
};

// --- Simulation State ---
let glucose = new Map(); let galactose = new Map(); let liveCells = new Map();
let frontierSet = new Set(); let start = false; let simulationSteps = 0; let lastStepBirths = 1;
let showGradientOverlay = false;

// --- Graphics & Grid ---
let nutrientBuffer, finalStateBuffer; let finalStateRendered = false; const SQRT3 = Math.sqrt(3); let hexWidth, hexHeight, vertSpacing;
const axialDirections = [ { q: +1, r: 0 }, { q: +1, r: -1 }, { q: 0, r: -1 }, { q: -1, r: 0 }, { q: -1, r: +1 }, { q: 0, r: +1 } ];
let canvasWidth, canvasHeight, bufferWidth, bufferHeight, gridOriginX, gridOriginY;

// --- Color Definitions ---
let nutrientColorG, nutrientColorH, nutrientColorMix; let cellColorH, cellColorG;

// --- UI Elements ---
let uiElements = { controlGroups: {} };
let uiParamDisplays = {};

// =============================================================================
// --- Math & Grid Helpers --- (No Changes)
// =============================================================================
function smoothstep(edge0, edge1, x) { x = constrain((x - edge0) / (edge1 - edge0), 0.0, 1.0); return x * x * (3 - 2 * x); }
function gaussianFalloff(dist, width) { if (width < epsilon) return (dist < epsilon) ? 1.0 : 0.0; return exp(-(dist * dist) / (2 * width * width)); }
function axialToCube(q, r) { return { x: q, y: -q - r, z: r }; }
function cubeToAxial(x, y, z) { return { q: x, r: z }; }
function cubeDistance(ax, ay, az, bx, by, bz) { return (Math.abs(ax - bx) + Math.abs(ay - by) + Math.abs(az - bz)) / 2; }
function axialDistance(q1, r1, q2, r2) { const ac = axialToCube(q1, r1); const bc = axialToCube(q2, r2); return cubeDistance(ac.x, ac.y, ac.z, bc.x, bc.y, bc.z); }
function hexToPixel(q, r) { let x = hexSize * (SQRT3 * q + SQRT3 / 2 * r); let y = hexSize * (3 / 2 * r); return { x: x, y: y }; }
function pixelToHexFractional(x, y) { let q = (SQRT3 / 3 * x - 1 / 3 * y) / hexSize; let r = (2 / 3 * y) / hexSize; return { q: q, r: r }; }
function cubeRound(x, y, z) { let rx = Math.round(x); let ry = Math.round(y); let rz = Math.round(z); let x_diff = Math.abs(rx - x); let y_diff = Math.abs(ry - y); let z_diff = Math.abs(rz - z); if (x_diff > y_diff && x_diff > z_diff) rx = -ry - rz; else if (y_diff > z_diff) ry = -rx - rz; else rz = -rx - ry; return { x: rx, y: ry, z: rz }; }
function pixelToHex(x, y) { const fracAxial = pixelToHexFractional(x, y); const fracCube = axialToCube(fracAxial.q, fracAxial.r); const roundedCube = cubeRound(fracCube.x, fracCube.y, fracCube.z); return cubeToAxial(roundedCube.x, roundedCube.y, roundedCube.z); }

// =============================================================================
// --- Drawing Helpers --- (No Changes)
// =============================================================================
function calculateGridDimensions() { /* ... */ }
function drawHexagon(buffer, q, r, fillColor, strokeColor = color(200, 180), weight = 0.5) { /* ... */ }

// =============================================================================
// --- P5.js Core Functions --- (No Changes)
// =============================================================================
function setup() {
    calculateGridDimensions();
    nutrientBuffer = createGraphics(bufferWidth, bufferHeight); finalStateBuffer = createGraphics(bufferWidth, bufferHeight); finalStateBuffer.pixelDensity(1);
    console.log(`Grid/Buffer Size: ${bufferWidth}x${bufferHeight}. Hex Size: ${hexSize}`);
    nutrientColorG = color(100, 220, 100, 200); nutrientColorH = color(200, 100, 220, 200); nutrientColorMix = color(240); cellColorH = color(60, 150, 240); cellColorG = color(255, 140, 40);
    let totalCanvasHeight = Math.max(canvasHeight, 900);
    let mainCanvas = createCanvas(canvasWidth, totalCanvasHeight); // Canvas width does NOT include UI panel
    let canvasContainer = select('#canvas-container'); if (canvasContainer) mainCanvas.parent('canvas-container'); else console.error("#canvas-container not found.");
    noiseDetail(4, 0.5);
    createUI();
    initGrid();
    background(240); push(); translate(gridOriginX, gridOriginY); image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2); drawLiveCells(this); pop();
    updateUIParamDisplays(); updateDynamicUIText(); noLoop();
}

function draw() {
     if (start) { stepSimulation(); }
    if (frameCount % drawInterval === 0 || !start || showGradientOverlay) {
        clear(); push(); translate(gridOriginX, gridOriginY);
        image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2);
        if (!showGradientOverlay) { if (start) { drawLiveCells(this); } else { if (finalStateRendered) image(finalStateBuffer, -bufferWidth / 2, -bufferHeight / 2); else drawLiveCells(this); } }
        pop();
    }
    updateDynamicUIText(); updateUIParamDisplays();
}

// =============================================================================
// --- Initialization & Reset ---
// =============================================================================
// Helper to calculate base G/H values for a given mode and position
function calculateBaseGradientValue(mode, q, r, R, pixelPos, minPixelX, maxPixelX, minPixelY, maxPixelY) {
    const { gradientCenter, gradientTransitionWidth,
            gaussianSourceQ, gaussianSourceR, gaussianWidth, gaussianPeakG, gaussianPeakH, gaussianBaseG, gaussianBaseH,
            stripeWidth // *** Added stripeWidth ***
          } = PARAMS;

    let baseG = 0.5, baseH = 0.5; // Default uniform mix

    if (mode === "smooth_radial") {
        let d = axialDistance(q, r, 0, 0); let maxDist = R > 0 ? R * 1.01 : 1; let d_norm = d / maxDist;
        let edge0 = gradientCenter - gradientTransitionWidth / 2; let edge1 = gradientCenter + gradientTransitionWidth / 2;
        let transitionFactor = smoothstep(edge0, edge1, d_norm); baseG = 1.0 - transitionFactor; baseH = transitionFactor;
    } else if (mode === "smooth_linear_h") {
        let d_norm = map(pixelPos.x, minPixelX, maxPixelX, 0, 1, true);
        let edge0 = gradientCenter - gradientTransitionWidth / 2; let edge1 = gradientCenter + gradientTransitionWidth / 2;
        let transitionFactor = smoothstep(edge0, edge1, d_norm); baseG = 1.0 - transitionFactor; baseH = transitionFactor;
    } else if (mode === "smooth_linear_v") {
        let d_norm = map(pixelPos.y, minPixelY, maxPixelY, 0, 1, true);
        let edge0 = gradientCenter - gradientTransitionWidth / 2; let edge1 = gradientCenter + gradientTransitionWidth / 2;
        let transitionFactor = smoothstep(edge0, edge1, d_norm); baseG = 1.0 - transitionFactor; baseH = transitionFactor;
    } else if (mode === "gaussian_source") {
        let d = axialDistance(q, r, gaussianSourceQ, gaussianSourceR); let distFactor = gaussianFalloff(d, gaussianWidth);
        baseG = gaussianBaseG + (gaussianPeakG - gaussianBaseG) * distFactor; baseH = gaussianBaseH + (gaussianPeakH - gaussianBaseH) * distFactor;
    }
    // *** NEW: Stripe Modes ***
      else if (mode === "stripes_h") {
          // Horizontal stripes based on pixel Y coordinate relative to buffer center
          let stripeIndex = floor((pixelPos.y + bufferHeight / 2) / stripeWidth) % 2;
          if (stripeIndex === 0) { baseG = 1.0; baseH = 0.0; } else { baseG = 0.0; baseH = 1.0; }
      } else if (mode === "stripes_v") {
          // Vertical stripes based on pixel X coordinate relative to buffer center
          let stripeIndex = floor((pixelPos.x + bufferWidth / 2) / stripeWidth) % 2;
          if (stripeIndex === 0) { baseG = 1.0; baseH = 0.0; } else { baseG = 0.0; baseH = 1.0; }
      }
    // *** End New Stripes ***
      else if (mode === "uniform_g") { baseG = 1.0; baseH = 0.0; }
      else if (mode === "uniform_h") { baseG = 0.0; baseH = 1.0; }
      // Default is uniform_mix (already set)

    return { g: constrain(baseG, 0, 1), h: constrain(baseH, 0, 1) };
}


function initGrid() {
    start = false; loop(); lastStepBirths = 1; finalStateRendered = false; showGradientOverlay = false;
    if(uiElements.gradientOverlayBtn) uiElements.gradientOverlayBtn.html('Show Gradient');
    if(finalStateBuffer) finalStateBuffer.clear(); simulationSteps = 0;
    glucose.clear(); galactose.clear(); liveCells.clear(); frontierSet.clear();

    const R = maxInitRadius; const bufferCenterX = bufferWidth / 2; const bufferCenterY = bufferHeight / 2;
    let minPixelX = Infinity, maxPixelX = -Infinity, minPixelY = Infinity, maxPixelY = -Infinity;
    for (let q = -R; q <= R; q++) { let r1 = Math.max(-R, -q - R); let r2 = Math.min(R, -q + R); for (let r = r1; r <= r2; r++) { const p = hexToPixel(q, r); minPixelX = Math.min(minPixelX, p.x); maxPixelX = Math.max(maxPixelX, p.x); minPixelY = Math.min(minPixelY, p.y); maxPixelY = Math.max(maxPixelY, p.y); } }
    const pixelWidth = (maxPixelX - minPixelX) || 1; const pixelHeight = (maxPixelY - minPixelY) || 1;

    nutrientBuffer.push(); nutrientBuffer.translate(bufferCenterX, bufferCenterY); nutrientBuffer.background(nutrientColorMix);
    let bufferNutrientG = color(red(nutrientColorG), green(nutrientColorG), blue(nutrientColorG));
    let bufferNutrientH = color(red(nutrientColorH), green(nutrientColorH), blue(nutrientColorH));

    for (let q = -R; q <= R; q++) {
        let r1 = Math.max(-R, -q - R); let r2 = Math.min(R, -q + R);
        for (let r = r1; r <= r2; r++) {
            const coordStr = `${q},${r}`; const pixelPos = hexToPixel(q, r); // Get pixel pos relative to grid 0,0
            const { gradientMode, noiseScale, noiseIntensityG, noiseIntensityH, noiseBaseGradient } = PARAMS;
            let gVal = 0.5, hVal = 0.5;

            // Calculate base value using helper function
            // Pass pixelPos relative to buffer center for stripe calcs if needed
            let baseValues = calculateBaseGradientValue(
                gradientMode === 'noisy' ? noiseBaseGradient : gradientMode,
                q, r, R, pixelPos, minPixelX, maxPixelX, minPixelY, maxPixelY
            );
            gVal = baseValues.g; hVal = baseValues.h;

            // Apply noise if selected
            if (gradientMode === 'noisy') {
                let noiseValG = noise(q * noiseScale, r * noiseScale, 0.1);
                let noiseValH = noise(q * noiseScale + 100, r * noiseScale + 100, 0.2);
                let noiseEffectG = (noiseValG - 0.5) * 2 * noiseIntensityG;
                let noiseEffectH = (noiseValH - 0.5) * 2 * noiseIntensityH;
                gVal += noiseEffectG; hVal += noiseEffectH;
            }

            gVal = constrain(gVal, 0, 1); hVal = constrain(hVal, 0, 1);
            glucose.set(coordStr, gVal); galactose.set(coordStr, hVal);

            let totalNutrient = gVal + hVal; let hexColor;
            if (totalNutrient < epsilon) hexColor = nutrientColorMix;
            else hexColor = lerpColor(bufferNutrientG, bufferNutrientH, hVal / totalNutrient);
            drawHexagon(nutrientBuffer, q, r, hexColor, null);
        }
    }
    nutrientBuffer.pop();
    console.log(`Grid Initialized. Grad: ${PARAMS.gradientMode}. Frontier Cleared.`);
    noLoop();
}


// =============================================================================
// --- Simulation Step Logic --- (Updated Galactose Fitness Calculation)
// =============================================================================
function stepSimulation() {
    // --- STOP CONDITION ---
    if (liveCells.size === 0 || (lastStepBirths === 0 && simulationSteps > 0)) { if (start) { /* ... Stop logic ... */ start = false; noLoop(); } return; }
    simulationSteps++;
    // *** Removed galactoseEfficiencyFactor from destructuring ***
    const { K_G, h_max, K_H, mutationProbability, minColonizationFitness, colonizationPressureFactor, maxColonizationProbability } = PARAMS;

    let nextLiveCells = new Map(liveCells); let nextFrontierSet = new Set(); let newbornCoords = []; let successfulBirthsCount = 0;

    // --- PASS 1: Calculate Fitness for all Current Frontier Cells ---
    let frontierFitness = new Map();
    for (const coordStr of frontierSet) {
        if (!liveCells.has(coordStr)) continue;
        const cell = liveCells.get(coordStr);
        const G = glucose.get(coordStr) || 0; const H = galactose.get(coordStr) || 0; let S = 0;
        const p_threshold = cell.p; const ratio = (G > epsilon) ? H / G : Infinity;
        const prefersGalactose = (ratio >= p_threshold && H > epsilon) || (G <= epsilon && H > epsilon);

        if (prefersGalactose) {
             // *** Simplified Galactose Fitness Calculation ***
             let s_potential_gal = h_max * (H / (K_H + H));
             // Scale potential by (1.0 - p), representing reduced efficiency for Glc specialists
             S = s_potential_gal * (1.0 - p_threshold);
        } else if (G > epsilon) {
            S = G / (K_G + G); // Glucose fitness remains the same
        }
        S = Math.max(0, S);
        frontierFitness.set(coordStr, S);

        const [q, r] = coordStr.split(',').map(Number); let hasEmpty = false;
        for(const dir of axialDirections){ const nq = q + dir.q; const nr = r + dir.r; const nCoordStr = `${nq},${nr}`; if(glucose.has(nCoordStr) && !liveCells.has(nCoordStr)){ hasEmpty = true; break; } }
        if(hasEmpty) nextFrontierSet.add(coordStr);
    }

    // --- PASS 2: Calculate Colonization Pressure & Chance for Each Relevant Empty Spot ---
    let potentialColonizations = new Map(); let emptySpotsChecked = new Set();
    for (const frontierCoordStr of frontierSet) { if (!frontierFitness.has(frontierCoordStr)) continue; const [q, r] = frontierCoordStr.split(',').map(Number); for (const dir of axialDirections) { const nq = q + dir.q; const nr = r + dir.r; const emptySpotCoord = `${nq},${nr}`; if (glucose.has(emptySpotCoord) && !liveCells.has(emptySpotCoord) && !emptySpotsChecked.has(emptySpotCoord)) { emptySpotsChecked.add(emptySpotCoord); let viableNeighbors = []; let totalPressure = 0; const [esq, esr] = emptySpotCoord.split(',').map(Number); for (const dir2 of axialDirections) { const nnq = esq + dir2.q; const nnr = esr + dir2.r; const neighborCoord = `${nnq},${nnr}`; if (liveCells.has(neighborCoord) && frontierFitness.has(neighborCoord)) { const neighborFitness = frontierFitness.get(neighborCoord); if (neighborFitness >= minColonizationFitness) { let pressure = map(neighborFitness, minColonizationFitness, 1.0, 0.01, 1.0, true); totalPressure += pressure; viableNeighbors.push({ coord: neighborCoord, fitness: neighborFitness, p: liveCells.get(neighborCoord).p }); } } } if (viableNeighbors.length > 0) { let P_spot_colonized = constrain(totalPressure * colonizationPressureFactor, 0, maxColonizationProbability); potentialColonizations.set(emptySpotCoord, { P_colonize: P_spot_colonized, viableNeighbors: viableNeighbors }); } } } }

    // --- PASS 3: Determine Actual Colonizations (Stochastic Check & Parent Selection) ---
    let colonizationEvents = [];
    for (const [spotCoord, data] of potentialColonizations.entries()) { if (random() < data.P_colonize) { let viableNeighbors = data.viableNeighbors; let totalCompetitorFitness = 0; for (const neighbor of viableNeighbors) totalCompetitorFitness += neighbor.fitness; if (totalCompetitorFitness > epsilon) { let randomThreshold = random(totalCompetitorFitness); let cumulativeFitness = 0; let winningParentData = viableNeighbors[viableNeighbors.length - 1]; for (const neighbor of viableNeighbors) { cumulativeFitness += neighbor.fitness; if (randomThreshold <= cumulativeFitness) { winningParentData = neighbor; break; } } colonizationEvents.push({ spotCoord: spotCoord, parentData: winningParentData }); } } }

    // --- PASS 4: Resolve Parent Conflicts (One offspring per parent) & Confirm Births ---
    shuffle(colonizationEvents, true); let parentsUsedThisStep = new Set();
    for (const event of colonizationEvents) { const spotCoord = event.spotCoord; const parentCoord = event.parentData.coord; if (!parentsUsedThisStep.has(parentCoord)) { let finalP = event.parentData.p; if (random() < mutationProbability) finalP = random(); const newCell = { p: finalP }; nextLiveCells.set(spotCoord, newCell); newbornCoords.push(spotCoord); successfulBirthsCount++; parentsUsedThisStep.add(parentCoord); } }

    // --- PASS 5: Update Frontier Set based on births ---
    for (const newbornCoordStr of newbornCoords) { nextFrontierSet.add(newbornCoordStr); const [q, r] = newbornCoordStr.split(',').map(Number); for (const dir of axialDirections) { const nq = q + dir.q; const nr = r + dir.r; const neighborCoordStr = `${nq},${nr}`; if (nextLiveCells.has(neighborCoordStr) && neighborCoordStr !== newbornCoordStr) { nextFrontierSet.add(neighborCoordStr); } } } nextFrontierSet = new Set([...nextFrontierSet].filter(coord => nextLiveCells.has(coord)));

    // --- PASS 6: Final State Update ---
    liveCells = nextLiveCells; frontierSet = nextFrontierSet; lastStepBirths = successfulBirthsCount; if (liveCells.size === 0) lastStepBirths = 0;
}

// =============================================================================
// --- Drawing Functions --- (No Changes)
// =============================================================================
function drawLiveCells(targetBuffer){ /* ... */ }
function drawLiveCellsOnBuffer(targetBuffer){ /* ... */ }

// =============================================================================
// --- UI Setup & Updates ---
// =============================================================================
function createUI() {
    // --- Helper Functions ---
    const createSliderRow = (label, paramKey, minVal, maxVal, stepVal, formatDigits, controlGroupArray) => {
        let rowDiv = createDiv('').parent('#ui-panel'); // Row container
        createDiv(label).html(label).parent(rowDiv); // Label
        let slider = createSlider(minVal, maxVal, PARAMS[paramKey], stepVal).parent(rowDiv)
             .input(() => { PARAMS[paramKey] = slider.value(); select('#' + paramKey + 'ValueDisplay').html(nf(slider.value(), 0, formatDigits)); if (paramKey.includes('gaussian') || paramKey.includes('gradient') || paramKey.includes('noise') || paramKey.includes('stripe')) { if (!start) { console.log(`Gradient param ${paramKey} changed, resetting...`); initGrid(); redraw(); updateUIParamDisplays(); } else { console.log(`Gradient param ${paramKey} changed, but sim running. Reset needed to apply.`); } } });
        createSpan(nf(PARAMS[paramKey], 1, formatDigits)).parent(rowDiv).id(paramKey + 'ValueDisplay');
        uiElements[paramKey + 'Slider'] = slider; uiParamDisplays[paramKey] = select('#' + paramKey + 'ValueDisplay');
        if (controlGroupArray) { controlGroupArray.push(rowDiv); } return rowDiv; };
    const createSectionHeader = (text, controlGroupArray) => { let headerP = createP(`<b>${text}</b>`).parent('#ui-panel'); if (controlGroupArray) { controlGroupArray.push(headerP); } return headerP; };

    // --- Create UI Elements ---
    createSectionHeader("Environment:");
    let gradTypeRowDiv = createDiv('').parent('#ui-panel').id('gradient-type-row');
    createSpan("Gradient Type:").parent(gradTypeRowDiv);
    uiElements.gradientSelector = createSelect().parent(gradTypeRowDiv)
        .changed(() => { PARAMS.gradientMode = uiElements.gradientSelector.value(); console.log("Gradient changed..."); if (!start) {initGrid(); redraw(); updateUIParamDisplays();} });
    // *** Added stripes modes ***
    ['smooth_radial', 'smooth_linear_h', 'smooth_linear_v', 'gaussian_source', 'noisy', 'stripes_h', 'stripes_v', 'uniform_g', 'uniform_h', 'uniform_mix'].forEach(opt => uiElements.gradientSelector.option(opt));
    uiElements.gradientSelector.selected(PARAMS.gradientMode);

    // -- Grouping Controls --
    uiElements.controlGroups.smooth_radial = []; uiElements.controlGroups.smooth_linear = []; uiElements.controlGroups.gaussian_source = []; uiElements.controlGroups.noisy = [];
    uiElements.controlGroups.stripes = []; // *** New group for stripes ***

    // --- Smooth Radial/Linear Controls ---
    createSectionHeader("Smooth Radial/Linear Controls:", uiElements.controlGroups.smooth_radial.concat(uiElements.controlGroups.smooth_linear));
    createSliderRow("Center (%):", 'gradientCenter', 0.0, 1.0, 0.01, 2, uiElements.controlGroups.smooth_radial.concat(uiElements.controlGroups.smooth_linear));
    createSliderRow("Transition Width (%):", 'gradientTransitionWidth', 0.05, 1.0, 0.01, 2, uiElements.controlGroups.smooth_radial.concat(uiElements.controlGroups.smooth_linear));

    // --- Gaussian Controls ---
    createSectionHeader("Gaussian Source Controls:", uiElements.controlGroups.gaussian_source);
    createSliderRow("Source Q:", 'gaussianSourceQ', -maxInitRadius, maxInitRadius, 1, 0, uiElements.controlGroups.gaussian_source);
    createSliderRow("Source R:", 'gaussianSourceR', -maxInitRadius, maxInitRadius, 1, 0, uiElements.controlGroups.gaussian_source);
    createSliderRow("Width (Hex):", 'gaussianWidth', 1, maxInitRadius/2, 0.5, 1, uiElements.controlGroups.gaussian_source);
    createSliderRow("Peak G:", 'gaussianPeakG', 0.0, 1.0, 0.01, 2, uiElements.controlGroups.gaussian_source);
    createSliderRow("Peak H:", 'gaussianPeakH', 0.0, 1.0, 0.01, 2, uiElements.controlGroups.gaussian_source);
    createSliderRow("Base G:", 'gaussianBaseG', 0.0, 1.0, 0.01, 2, uiElements.controlGroups.gaussian_source);
    createSliderRow("Base H:", 'gaussianBaseH', 0.0, 1.0, 0.01, 2, uiElements.controlGroups.gaussian_source);

    // --- Noisy Controls ---
    createSectionHeader("Noisy Gradient Controls:", uiElements.controlGroups.noisy);
    createSliderRow("Noise Scale:", 'noiseScale', 0.01, 0.5, 0.01, 2, uiElements.controlGroups.noisy);
    createSliderRow("Noise Inten. G:", 'noiseIntensityG', 0.0, 0.5, 0.01, 2, uiElements.controlGroups.noisy);
    createSliderRow("Noise Inten. H:", 'noiseIntensityH', 0.0, 0.5, 0.01, 2, uiElements.controlGroups.noisy);
    let noiseBaseRowDiv = createDiv('').parent('#ui-panel').id('noise-base-row');
    createSpan("Noise Base Grad:").parent(noiseBaseRowDiv);
    uiElements.noiseBaseSelector = createSelect().parent(noiseBaseRowDiv) .changed(() => { PARAMS.noiseBaseGradient = uiElements.noiseBaseSelector.value(); if (!start) { initGrid(); redraw(); } });
    ['smooth_radial', 'smooth_linear_h', 'smooth_linear_v', 'gaussian_source', 'uniform_g', 'uniform_h', 'uniform_mix'].forEach(opt => uiElements.noiseBaseSelector.option(opt));
    uiElements.noiseBaseSelector.selected(PARAMS.noiseBaseGradient);
    uiElements.controlGroups.noisy.push(noiseBaseRowDiv);

    // --- *** Stripe Controls *** ---
    createSectionHeader("Stripe Controls:", uiElements.controlGroups.stripes);
    createSliderRow("Stripe Width (px):", 'stripeWidth', 5, 100, 1, 0, uiElements.controlGroups.stripes);


    // --- Cell Parameters Section ---
    createSectionHeader("Cell Parameters:");
    createSliderRow("Mutation Prob:", 'mutationProbability', 0, 0.1, 0.001, 3);
    createSliderRow("Initial P (H/G Thr):", 'initialP', 0, 1, 0.01, 2);
    createSliderRow("Glc Affinity (K<sub>G</sub>):", 'K_G', 0.01, 1.0, 0.01, 2);
    createSliderRow("Gal Affinity (K<sub>H</sub>):", 'K_H', 0.01, 1.0, 0.01, 2);
    // *** Updated Label for h_max ***
    createSliderRow("MAX Gal Rate (h<sub>max</sub>):", 'h_max', 0.1, 1.5, 0.01, 2);
    // *** REMOVED Gal Max Efficiency Slider ***

    // --- Colonization Parameters Section ---
    createSectionHeader("Colonization Parameters:");
    createSliderRow("Min Viability Fitness:", 'minColonizationFitness', 0, 0.5, 0.01, 2);
    createSliderRow("Pressure Factor:", 'colonizationPressureFactor', 0.1, 2.0, 0.1, 1);
    createSliderRow("Max Colonize Prob:", 'maxColonizationProbability', 0.5, 1.0, 0.01, 2);

    // --- Controls Section ---
    createSectionHeader("Controls:");
    let controlsDiv = createDiv('').parent('#ui-panel').id('controls-container');
    uiElements.startBtn = createButton("Start / Pause (Space)").parent(controlsDiv).mousePressed(toggleSimulation);
    uiElements.resetBtn = createButton("Reset Grid").parent(controlsDiv).mousePressed(() => { console.log("Resetting grid..."); initGrid(); redraw(); updateDynamicUIText(); });
    uiElements.gradientOverlayBtn = createButton("Show Gradient").parent(controlsDiv).mousePressed(toggleGradientOverlay);

    // --- Status Display Area ---
    let statusDiv = createDiv('').parent('#ui-panel').id('status-section');
    createP("<b>Status:</b>").parent(statusDiv);
    createSpan('Step: 0').parent(statusDiv).id('stepValueDisplay');
    createSpan('Cell Count: 0').parent(statusDiv).id('cellCountValueDisplay');
    createSpan('Status: Initialized').parent(statusDiv).id('statusValueDisplay');

    updateUIParamDisplays(); // Initial setup of visibility
}

function toggleSimulation() { /* ... (no change) ... */ }
function toggleGradientOverlay() { /* ... (no change) ... */ }

function updateUIParamDisplays() {
    // Update slider values displayed next to them
    for (const key in uiParamDisplays) { /* ... (Update logic remains same) ... */ }

    // Update visibility of gradient controls
    const currentMode = PARAMS.gradientMode;
    for (const groupName in uiElements.controlGroups) {
        let showGroup = false;
        if (groupName === 'smooth_radial' && currentMode === 'smooth_radial') showGroup = true;
        else if (groupName === 'smooth_linear' && (currentMode === 'smooth_linear_h' || currentMode === 'smooth_linear_v')) showGroup = true;
        else if (groupName === 'gaussian_source' && currentMode === 'gaussian_source') showGroup = true;
        else if (groupName === 'noisy' && currentMode === 'noisy') showGroup = true;
        // *** Added check for stripes ***
        else if (groupName === 'stripes' && (currentMode === 'stripes_h' || currentMode === 'stripes_v')) showGroup = true;

        uiElements.controlGroups[groupName].forEach(element => { if (element) { if (showGroup) element.show(); else element.hide(); } });
    }
}

function updateDynamicUIText() { /* ... (no change) ... */ }

// =============================================================================
// --- User Interaction --- (No Changes Needed)
// =============================================================================
function mousePressed() { /* ... */ }
function plantPatch(centerQ, centerR, radius) { /* ... */ }
function keyPressed() { /* ... */ }
// =============================================================================
