// --- Core Microbial Strategy Simulation (v31 - Julia Model Port Refined: Random Alternating Bands & Simplified UI) ---
// Features:
// - Nutrient Environment: RandomHexBandNutrientProfile (alternating Glu-High/Gal-High bands)
//   - Dynamic generation of bands as colony expands.
//   - Initial "homeland" environment.
// - Cell Model: Genotype (Red/Green), Galactose Phenotype (ON/OFF), Lag phase with memory.
// - Simplified UI: Fewer parameters exposed, more defaults.
// - Colonization Model (p5.js style): Viability Threshold, Collective Pressure, Probabilistic Success, Parent Limitation.
// - User interaction: Click to plant, Spacebar Start/Pause, Gradient Overlay Toggle.
// --- NO CELL DEATH in this version ---
// ---

// --- Simulation Settings ---
let maxWorldRadius = 60;
let hexSize = 5;
const epsilon = 1e-9;
let drawInterval = 1; // Faster drawing for dynamic env

// --- MODEL PARAMETERS (Defaults, adjustable via UI where exposed) ---
let PARAMS = {
  // --- Nutrient Environment: Random Alternating Bands (NEW) ---
  nutrientProfile: {
    R0_homeland_radius: 1.0,
    homeland_env_type: "GluHigh_GalLow", // "GluHigh_GalLow" or "GluLow_GalHigh"
    min_band_width: 2.0,
    max_band_width: 7.0,
    // Specific nutrient levels for the two environment types:
    env_types: {
      GluHigh_GalLow: { glu: 0.9, gal: 0.1 },
      GluLow_GalHigh: { glu: 0.1, gal: 0.9 }
    },
    generation_lookahead_radius: 10.0 // How far ahead to generate bands
  },

  // --- Cell Genotype, Phenotype & Lag (Simplified exposure) ---
  initialGenotypeMix: 0.5, // Fraction of Green (1) cells if initial_condition_type is 'mix'
  initialConditionType: "mix", // "mix", "all_green", "all_red"
  mutationRateGreenToRed: 0.0005,
  mutationRateRedToGreen: 0.00005,

  maxLagTime: 0.6,
  memoryDecayGens: 4.0,
  minLagTimeWithMemory: 0.05,
  lagTimeStochasticity: 0.1, // StdDev fraction for lag time noise (0 for no noise)
  dtPerCell: 0.02,     // Effective time step advanced for lag reduction per cell added

  // --- Fitness Parameters (Mostly fixed, based on Julia example) ---
  fitnessParams: {
    fitness_green_on_gal: 1.3,
    fitness_red_on_gal: 0.7,
    fitness_on_no_gal_cost_factor: 0.8,
    fitness_green_off_glu: 1.0,
    fitness_red_off_glu: 1.0, // Red equally good on glucose
    fitness_off_on_gal: 0.1,
    fitness_during_lag: 0.005,
  },
  // Nutrient thresholds for phenotype switching & fitness
  nutrientThresholdForPresence: 0.4, // Min nutrient level to be "present"
  glucoseThresholdForRepression: 0.6, // High Glc level that represses Gal

  // --- Colonization ---
  minColonizationFitness: 0.05,
  colonizationPressureFactor: 0.8,
  maxColonizationProbability: 0.95,
  initialPatchRadius: 1, // Radius of patch planted on click (smaller default)

  // --- OLD Gradient Modes (accessible if gradientMode changed from UI, but not default) ---
  // Defaulting to 'random_alternating_bands'
  gradientMode: "random_alternating_bands",
  // Smooth Radial/Linear Controls
  gradientCenter: 0.25, gradientTransitionWidth: 0.5,
  // Gaussian Source Controls
  gaussianSourceQ: 0, gaussianSourceR: 0, gaussianWidth: 15,
  gaussianPeakG: 1.0, gaussianPeakH: 0.0, gaussianBaseG: 0.0, gaussianBaseH: 0.2,
  // Noisy Gradient Controls
  noiseScale: 0.1, noiseIntensityG: 0.3, noiseIntensityH: 0.3, noiseBaseGradient: "uniform_mix",
  // Stripe Controls
  stripeWidth: 40,
};

// --- Simulation State ---
let glucose = new Map(); let galactose = new Map(); let liveCells = new Map();
let frontierSet = new Set(); let start = false; let simulationSteps = 0; let lastStepBirths = 1;
let showGradientOverlay = false;

// --- Nutrient Profile State (NEW) ---
// Structure: [{r_in: val, r_out: val, type: "GluHigh_GalLow" / "GluLow_GalHigh"}, ...]
let currentNutrientBands = [];
let currentMaxGeneratedRadius_Nutrients = 0;
let currentColonyMaxRadius = 0; // Actual max radius of living cells

// --- Graphics & Grid ---
let nutrientBuffer, finalStateBuffer; let finalStateRendered = false; const SQRT3 = Math.sqrt(3); let hexWidth, hexHeight, vertSpacing;
const axialDirections = [ { q: +1, r: 0 }, { q: +1, r: -1 }, { q: 0, r: -1 }, { q: -1, r: 0 }, { q: -1, r: +1 }, { q: 0, r: +1 } ];
let canvasWidth, canvasHeight, bufferWidth, bufferHeight, gridOriginX, gridOriginY;

// --- Color Definitions ---
// Nutrient environment band colors (for overlay)
let envColorGluHighGalLow, envColorGluLowGalHigh, envColorMix;
// Cell Colors (Genotype & Phenotype) - Updated based on Julia plot
let cellColor_Green_GalOFF, cellColor_Green_GalON, cellColor_Green_Lag;
let cellColor_Red_GalOFF, cellColor_Red_GalON, cellColor_Red_Lag;

// --- UI Elements ---
let uiElements = { controlGroups: {} };
let uiParamDisplays = {};

// =============================================================================
// --- Math & Grid Helpers (Unchanged from v30) ---
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
// --- Drawing Helpers (Unchanged) ---
// =============================================================================
function calculateGridDimensions() {
    hexWidth = SQRT3 * hexSize; hexHeight = 2 * hexSize; vertSpacing = hexHeight * 3/4;
    bufferWidth = Math.ceil((maxWorldRadius * 2 + 5) * hexWidth); bufferHeight = Math.ceil((maxWorldRadius * 2 + 5) * vertSpacing);
    canvasWidth = bufferWidth + 40; canvasHeight = bufferHeight + 40;
    gridOriginX = canvasWidth / 2; gridOriginY = canvasHeight / 2;
}

function drawHexagon(buffer, q, r, fillColor, strokeColor = color(200, 180), weight = 0.5) {
    const p = hexToPixel(q, r); buffer.fill(fillColor);
    if (strokeColor) { buffer.stroke(strokeColor); buffer.strokeWeight(weight); }
    else { buffer.noStroke(); }
    buffer.beginShape();
    for (let i = 0; i < 6; i++) { let angle = Math.PI / 3 * i + Math.PI / 6; buffer.vertex(p.x + hexSize * Math.cos(angle), p.y + hexSize * Math.sin(angle)); }
    buffer.endShape(CLOSE);
 }

// =============================================================================
// --- P5.js Core Functions ---
// =============================================================================
function setup() {
    calculateGridDimensions();
    nutrientBuffer = createGraphics(bufferWidth, bufferHeight); finalStateBuffer = createGraphics(bufferWidth, bufferHeight); finalStateBuffer.pixelDensity(1);
    console.log(`Grid/Buffer Size: ${bufferWidth}x${bufferHeight}. Hex Size: ${hexSize}`);

    // Nutrient Environment Colors (for gradient display)
    envColorGluHighGalLow = color(128, 0, 128, 150); // Purpleish (mediumpurple in Julia)
    envColorGluLowGalHigh = color(0, 128, 128, 150); // Tealish (lightseagreen in Julia)
    envColorMix         = color(211, 211, 211, 150); // Lightgray (for unassigned or old modes)

    // Cell Colors (matching Julia plot)
    cellColor_Green_GalOFF = color(152, 251, 152); // palegreen
    cellColor_Green_GalON  = color(34, 139, 34);   // forestgreen
    cellColor_Green_Lag    = color(255, 215, 0);   // gold
    cellColor_Red_GalOFF   = color(255, 182, 193); // lightpink
    cellColor_Red_GalON    = color(128, 0, 0);     // maroon
    cellColor_Red_Lag      = color(255, 140, 0);   // darkorange

    let totalCanvasHeight = Math.max(canvasHeight, 750); // Adjusted for simpler UI
    let mainCanvas = createCanvas(canvasWidth, totalCanvasHeight);
    let canvasContainer = select('#canvas-container'); if (canvasContainer) mainCanvas.parent('canvas-container'); else console.error("#canvas-container not found.");
    noiseDetail(4, 0.5);
    createUI();
    initGridAndNutrients();
    background(240); push(); translate(gridOriginX, gridOriginY); image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2); drawLiveCells(this); pop();
    updateUIParamDisplays(); updateDynamicUIText(); noLoop();
}

function draw() {
     if (start) { stepSimulation(); }
    if (frameCount % drawInterval === 0 || !start || showGradientOverlay || (PARAMS.gradientMode === "random_alternating_bands" && start)) {
        clear(); push(); translate(gridOriginX, gridOriginY);
        image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2);
        if (!showGradientOverlay) {
             if (start) { drawLiveCells(this); }
             else {
                 if (finalStateRendered) image(finalStateBuffer, -bufferWidth / 2, -bufferHeight / 2);
                 else drawLiveCells(this);
             }
        }
        pop();
    }
    updateDynamicUIText(); // updateUIParamDisplays() only needed if PARAMS change by other means
}

// =============================================================================
// --- Nutrient Environment Logic (NEW - Random Alternating Bands) ---
// =============================================================================

function ensureNutrientBandsGenerated() {
    const np = PARAMS.nutrientProfile;
    if (currentColonyMaxRadius + np.generation_lookahead_radius > currentMaxGeneratedRadius_Nutrients) {
        let nextBandStartRadius = currentMaxGeneratedRadius_Nutrients;
        let lastEnvType;

        if (currentNutrientBands.length === 0) { // Initial generation or reset
            lastEnvType = (np.homeland_env_type === "GluHigh_GalLow") ? "GluLow_GalHigh" : "GluHigh_GalLow"; // So the first *new* band alternates
        } else {
            lastEnvType = currentNutrientBands[currentNutrientBands.length - 1].type;
        }

        let newBandsAdded = false;
        while (currentMaxGeneratedRadius_Nutrients < currentColonyMaxRadius + np.generation_lookahead_radius) {
            let currentEnvTypeToAdd = (lastEnvType === "GluHigh_GalLow") ? "GluLow_GalHigh" : "GluHigh_GalLow";
            let bandWidth = random(np.min_band_width, np.max_band_width);
            let bandEndRadius = nextBandStartRadius + bandWidth;

            // Safety: ensure bands don't exceed maxWorldRadius by too much if colony is huge
            bandEndRadius = min(bandEndRadius, maxWorldRadius * 1.2);
            if (bandEndRadius <= nextBandStartRadius + 0.1) break; // Avoid tiny/zero width bands near edge

            currentNutrientBands.push({
                r_in: nextBandStartRadius,
                r_out: bandEndRadius,
                type: currentEnvTypeToAdd
            });
            // console.log(`Nutrient band: R=[${nextBandStartRadius.toFixed(1)}, ${bandEndRadius.toFixed(1)}), Type: ${currentEnvTypeToAdd}`);
            newBandsAdded = true;

            currentMaxGeneratedRadius_Nutrients = bandEndRadius;
            nextBandStartRadius = bandEndRadius;
            lastEnvType = currentEnvTypeToAdd;

            if (currentMaxGeneratedRadius_Nutrients >= maxWorldRadius * 1.2) break; // Stop if we've covered well beyond world
        }
        if (newBandsAdded) {
             console.log(`Nutrient bands extended. Current Max Gen R: ${currentMaxGeneratedRadius_Nutrients.toFixed(1)}. Colony Max R: ${currentColonyMaxRadius.toFixed(1)}`);
        }
        return newBandsAdded; // Return true if bands were actually added/changed
    }
    return false; // No change to bands needed
}


function getNutrientLevelsAtHex(q, r) {
    const np = PARAMS.nutrientProfile; // Current nutrient profile settings
    const dist = axialDistance(q, r, 0, 0);

    if (PARAMS.gradientMode === "random_alternating_bands") {
        // Check homeland first
        if (dist <= np.R0_homeland_radius) {
            return np.env_types[np.homeland_env_type];
        }
        // Check generated bands
        for (const band of currentNutrientBands) {
            if (dist > band.r_in && dist <= band.r_out) { // Note: > r_in to be exclusive of previous band's r_out
                return np.env_types[band.type];
            }
        }
        // If beyond generated bands (should be rare if ensureNutrientBandsGenerated is called)
        // Default to the type of the "next" band that *would* be generated
        let lastKnownType = np.homeland_env_type;
        if (currentNutrientBands.length > 0) {
            lastKnownType = currentNutrientBands[currentNutrientBands.length - 1].type;
        }
        const nextType = (lastKnownType === "GluHigh_GalLow") ? "GluLow_GalHigh" : "GluHigh_GalLow";
        // console.warn(`Radius ${dist.toFixed(1)} for hex (${q},${r}) outside defined nutrient bands. Defaulting to next type: ${nextType}. MaxGenR: ${currentMaxGeneratedRadius_Nutrients.toFixed(1)}`);
        return np.env_types[nextType];

    } else { // Fallback to old static gradient logic from v30
        const R_world = maxWorldRadius;
        const pixelPos = hexToPixel(q, r);
        const approxMinPixelX = -R_world * hexWidth / 2;
        const approxMaxPixelX = R_world * hexWidth / 2;
        const approxMinPixelY = -R_world * vertSpacing / 2;
        const approxMaxPixelY = R_world * vertSpacing / 2;

        let baseValues = calculateBaseGradientValue_Static(
            PARAMS.gradientMode === 'noisy' ? PARAMS.noiseBaseGradient : PARAMS.gradientMode,
            q, r, R_world, pixelPos, approxMinPixelX, approxMaxPixelX, approxMinPixelY, approxMaxPixelY
        );
        let gVal = baseValues.g; let hVal = baseValues.h;

        if (PARAMS.gradientMode === 'noisy') {
            let noiseValG = noise(q * PARAMS.noiseScale, r * PARAMS.noiseScale, 0.1);
            let noiseValH = noise(q * PARAMS.noiseScale + 100, r * PARAMS.noiseScale + 100, 0.2);
            gVal += (noiseValG - 0.5) * 2 * PARAMS.noiseIntensityG;
            hVal += (noiseValH - 0.5) * 2 * PARAMS.noiseIntensityH;
        }
        return { g: constrain(gVal, 0, 1), h: constrain(hVal, 0, 1) };
    }
}

// v30's static gradient calculator, kept for other modes
function calculateBaseGradientValue_Static(mode, q, r, R, pixelPos, minPixelX, maxPixelX, minPixelY, maxPixelY) {
    const { gradientCenter, gradientTransitionWidth,
            gaussianSourceQ, gaussianSourceR, gaussianWidth, gaussianPeakG, gaussianPeakH, gaussianBaseG, gaussianBaseH,
            stripeWidth
          } = PARAMS; // Uses global PARAMS
    let baseG = 0.5, baseH = 0.5;
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
    } else if (mode === "stripes_h") {
          let stripeIndex = floor((pixelPos.y + bufferHeight / 2) / stripeWidth) % 2;
          if (stripeIndex === 0) { baseG = 1.0; baseH = 0.0; } else { baseG = 0.0; baseH = 1.0; }
      } else if (mode === "stripes_v") {
          let stripeIndex = floor((pixelPos.x + bufferWidth / 2) / stripeWidth) % 2;
          if (stripeIndex === 0) { baseG = 1.0; baseH = 0.0; } else { baseG = 0.0; baseH = 1.0; }
      } else if (mode === "uniform_g") { baseG = 1.0; baseH = 0.0; }
      else if (mode === "uniform_h") { baseG = 0.0; baseH = 1.0; }
      else if (mode === "uniform_mix") { /* already 0.5, 0.5 */ }
    return { g: constrain(baseG, 0, 1), h: constrain(baseH, 0, 1) };
}


function populateNutrientMapsAndDrawBuffer() {
    const R_world = maxWorldRadius;
    const bufferCenterX = bufferWidth / 2;
    const bufferCenterY = bufferHeight / 2;

    nutrientBuffer.push();
    nutrientBuffer.translate(bufferCenterX, bufferCenterY);
    nutrientBuffer.background(envColorMix); // Default background

    for (let q_ax = -R_world; q_ax <= R_world; q_ax++) {
        let r1 = Math.max(-R_world, -q_ax - R_world);
        let r2 = Math.min(R_world, -q_ax + R_world);
        for (let r_ax = r1; r_ax <= r2; r_ax++) {
            const coordStr = `${q_ax},${r_ax}`;
            const { g, h } = getNutrientLevelsAtHex(q_ax, r_ax); // This now uses the band logic
            glucose.set(coordStr, g);
            galactose.set(coordStr, h);

            let hexColor = envColorMix; // Default
            if (PARAMS.gradientMode === "random_alternating_bands") {
                // Color based on the dominant nutrient type from the band definition
                const dist = axialDistance(q_ax, r_ax, 0, 0);
                let envType = PARAMS.nutrientProfile.homeland_env_type; // Start with homeland
                if (dist > PARAMS.nutrientProfile.R0_homeland_radius) {
                    let bandFound = false;
                    for (const band of currentNutrientBands) {
                        if (dist > band.r_in && dist <= band.r_out) {
                            envType = band.type;
                            bandFound = true;
                            break;
                        }
                    }
                    if (!bandFound) { // If outside, determine next type
                        let lastKnownType = PARAMS.nutrientProfile.homeland_env_type;
                        if (currentNutrientBands.length > 0) {
                           lastKnownType = currentNutrientBands[currentNutrientBands.length - 1].type;
                        }
                        envType = (lastKnownType === "GluHigh_GalLow") ? "GluLow_GalHigh" : "GluHigh_GalLow";
                    }
                }
                hexColor = (envType === "GluHigh_GalLow") ? envColorGluHighGalLow : envColorGluLowGalHigh;
            } else { // For old static modes, color based on G/H ratio
                if (g > h && g > PARAMS.nutrientThresholdForPresence * 0.5) hexColor = lerpColor(envColorMix, envColorGluHighGalLow, g);
                else if (h >= g && h > PARAMS.nutrientThresholdForPresence * 0.5) hexColor = lerpColor(envColorMix, envColorGluLowGalHigh, h);
            }
            drawHexagon(nutrientBuffer, q_ax, r_ax, hexColor, null);
        }
    }
    nutrientBuffer.pop();
    if (start) redraw();
}


// =============================================================================
// --- Initialization & Reset ---
// =============================================================================
function initGridAndNutrients() {
    start = false; loop(); lastStepBirths = 1; finalStateRendered = false; showGradientOverlay = false;
    if(uiElements.buttons?.gradientOverlayBtn) uiElements.buttons.gradientOverlayBtn.html('Show Gradient');
    if(finalStateBuffer) finalStateBuffer.clear(); simulationSteps = 0;
    glucose.clear(); galactose.clear(); liveCells.clear(); frontierSet.clear();

    currentColonyMaxRadius = PARAMS.nutrientProfile.R0_homeland_radius; // Initialize with homeland
    currentNutrientBands = []; // Clear existing bands
    currentMaxGeneratedRadius_Nutrients = PARAMS.nutrientProfile.R0_homeland_radius; // Start generation from homeland edge
    
    if (PARAMS.gradientMode === "random_alternating_bands") {
        ensureNutrientBandsGenerated(); // Generate initial bands beyond homeland
    }
    populateNutrientMapsAndDrawBuffer(); // Populates G/H maps and draws buffer

    console.log(`Grid Initialized. Mode: ${PARAMS.gradientMode}. Homeland R: ${PARAMS.nutrientProfile.R0_homeland_radius}, Type: ${PARAMS.nutrientProfile.homeland_env_type}`);
    noLoop();
    updateDynamicUIText();
}

// =============================================================================
// --- Cell Fitness Logic (Slightly adjusted thresholds from Julia) ---
// =============================================================================
function calculateCellFitness(cellData, G, H, fitnessP) { // fitnessP is PARAMS.fitnessParams
    if (cellData.lagTimeRemaining > epsilon) {
        return fitnessP.fitness_during_lag;
    }

    const { nutrientThresholdForPresence, glucoseThresholdForRepression } = PARAMS;
    let isGalPresent = H >= nutrientThresholdForPresence; // Use >= for thresholds
    let isGluHighForRepression = G >= glucoseThresholdForRepression;

    if (cellData.galPhenotype === "Gal_ON") {
        if (isGalPresent && !isGluHighForRepression) {
            return cellData.genotype === 1 ? fitnessP.fitness_green_on_gal : fitnessP.fitness_red_on_gal;
        } else { // Gal_ON, but Gal is low/absent OR Glu is high (repressive)
            // Assumes cell tries to use Glucose if Gal pathway is "on" but non-productive
            let base_glu_fitness = (G >= nutrientThresholdForPresence) ?
                                   (cellData.genotype === 1 ? fitnessP.fitness_green_off_glu : fitnessP.fitness_red_off_glu)
                                   : 0.001; // Very low if no alternative
            return base_glu_fitness * fitnessP.fitness_on_no_gal_cost_factor;
        }
    } else { // Gal_OFF
        if (isGalPresent && !isGluHighForRepression) { // Should be inducing, but currently OFF
            return fitnessP.fitness_off_on_gal;
        } else { // Gal_OFF, and Gal is low/absent OR Glu is present (normal Glu metabolism or low nutrient)
            return (G >= nutrientThresholdForPresence) ?
                   (cellData.genotype === 1 ? fitnessP.fitness_green_off_glu : fitnessP.fitness_red_off_glu)
                   : 0.001; // Very low if no Glucose
        }
    }
}


// =============================================================================
// --- Simulation Step Logic (Largely v30, adapted for new nutrient gen) ---
// =============================================================================
function stepSimulation() {
    if (liveCells.size === 0 || (lastStepBirths === 0 && simulationSteps > 0)) {
        if (start) {
            let reason = liveCells.size === 0 ? "Extinction" : "Stasis";
            console.log(`${reason} @ step ${simulationSteps}. Stopping.`);
            if (liveCells.size > 0 && !finalStateRendered) {
                drawLiveCellsOnBuffer(finalStateBuffer); finalStateRendered = true;
            }
            start = false; noLoop(); updateDynamicUIText();
        }
        return;
    }
    simulationSteps++;

    const { minColonizationFitness, colonizationPressureFactor, maxColonizationProbability,
            mutationRateGreenToRed, mutationRateRedToGreen,
            maxLagTime, memoryDecayGens, minLagTimeWithMemory, lagTimeStochasticity,
            dtPerCell, fitnessParams, nutrientThresholdForPresence, glucoseThresholdForRepression
          } = PARAMS;

    let nextLiveCells = new Map(liveCells);
    let nextFrontierSet = new Set();
    let newbornCoordsAndData = [];
    let successfulBirthsCount = 0;

    // --- Update currentColonyMaxRadius (based on ALL live cells, not just frontier) ---
    let tempMaxR = PARAMS.nutrientProfile.R0_homeland_radius; // Start with homeland min
    if (liveCells.size > 0) {
        for (const [coordStr, cell] of liveCells) {
            const [q,r] = coordStr.split(',').map(Number);
            tempMaxR = max(tempMaxR, axialDistance(q,r,0,0));
        }
    }
    currentColonyMaxRadius = tempMaxR;


    // --- PASS 1: Calculate Fitness for all Current Frontier Cells ---
    let frontierFitness = new Map();
    for (const coordStr of frontierSet) {
        if (!liveCells.has(coordStr)) continue;
        const cell = liveCells.get(coordStr);
        const G = glucose.get(coordStr) || 0;
        const H = galactose.get(coordStr) || 0;
        let fitness = calculateCellFitness(cell, G, H, fitnessParams);
        frontierFitness.set(coordStr, fitness);

        const [q, r] = coordStr.split(',').map(Number); let hasEmpty = false;
        for(const dir of axialDirections){ const nq = q + dir.q; const nr = r + dir.r; if(glucose.has(`${nq},${nr}`) && !liveCells.has(`${nq},${nr}`)){ hasEmpty = true; break; } }
        if(hasEmpty) nextFrontierSet.add(coordStr);
    }

    // --- PASS 2 & 3: Colonization Pressure, Stochastic Success (v30 logic) ---
    let potentialColonizations = new Map(); let emptySpotsChecked = new Set();
    for (const frontierCoordStr of frontierSet) {
        if (!frontierFitness.has(frontierCoordStr)) continue;
        const [q, r] = frontierCoordStr.split(',').map(Number);
        for (const dir of axialDirections) {
            const nq = q + dir.q; const nr = r + dir.r;
            const emptySpotCoord = `${nq},${nr}`;
            if (glucose.has(emptySpotCoord) && !liveCells.has(emptySpotCoord) && !emptySpotsChecked.has(emptySpotCoord)) {
                emptySpotsChecked.add(emptySpotCoord);
                let viableNeighbors = []; let totalPressure = 0;
                const [esq, esr] = emptySpotCoord.split(',').map(Number);
                for (const dir2 of axialDirections) {
                    const nnq = esq + dir2.q; const nnr = esr + dir2.r;
                    const neighborCoord = `${nnq},${nnr}`;
                    if (liveCells.has(neighborCoord) && frontierFitness.has(neighborCoord)) {
                        const neighborFitness = frontierFitness.get(neighborCoord);
                        if (neighborFitness >= minColonizationFitness) {
                            totalPressure += map(neighborFitness, minColonizationFitness, 1.0, 0.01, 1.0, true);
                            viableNeighbors.push({ coord: neighborCoord, fitness: neighborFitness, cellData: liveCells.get(neighborCoord) });
                        }
                    }
                }
                if (viableNeighbors.length > 0) {
                    potentialColonizations.set(emptySpotCoord, { P_colonize: constrain(totalPressure * colonizationPressureFactor, 0, maxColonizationProbability), viableNeighbors: viableNeighbors });
                }
            }
        }
    }
    let colonizationEvents = [];
    for (const [spotCoord, data] of potentialColonizations.entries()) { if (random() < data.P_colonize) colonizationEvents.push({ spotCoord: spotCoord, viableNeighbors: data.viableNeighbors }); }


    // --- PASS 4: Resolve Parent Conflicts, Determine Daughter Genotype/Phenotype, Confirm Births (v30 logic) ---
    shuffle(colonizationEvents, true);
    let parentsUsedThisStep = new Set();

    for (const event of colonizationEvents) {
        const spotCoord = event.spotCoord;
        const viableParents = event.viableNeighbors;
        let totalParentFitness = 0; viableParents.forEach(p => totalParentFitness += p.fitness);
        if (totalParentFitness < epsilon) continue;

        let chosenPhenoParentData = null, phenoParentCoord = null, foundPhenoParent = false;
        let shuffledViableParents = shuffle([...viableParents]);
        for (const parent of shuffledViableParents) { // Try to pick pheno parent
            if (!parentsUsedThisStep.has(parent.coord)) {
                // Simplified: first available parent in shuffled list is chosen if it has some fitness
                // A more complex weighted choice could be re-added if needed.
                if (parent.fitness > epsilon) {
                    chosenPhenoParentData = parent.cellData;
                    phenoParentCoord = parent.coord;
                    foundPhenoParent = true;
                    break;
                }
            }
        }
        if (!foundPhenoParent) continue; // All potential parents already reproduced or had zero fitness

        let weightedSumGreenContrib = 0.0, weightedSumRedContrib = 0.0;
        for (const p of viableParents) { if (p.cellData.genotype === 1) weightedSumGreenContrib += p.fitness; else weightedSumRedContrib += p.fitness; }
        let pG_genotype = (weightedSumGreenContrib + weightedSumRedContrib > epsilon) ? (weightedSumGreenContrib / (weightedSumGreenContrib + weightedSumRedContrib)) : (weightedSumGreenContrib > 0 ? 1.0 : 0.0);
        let daughterGenotype = random() < pG_genotype ? 1 : 0;
        if (daughterGenotype === 1 && random() < mutationRateGreenToRed) daughterGenotype = 0;
        if (daughterGenotype === 0 && random() < mutationRateRedToGreen) daughterGenotype = 1;

        let daughterCell = {
            genotype: daughterGenotype,
            galPhenotype: chosenPhenoParentData.galPhenotype,
            gensSinceGalOn: chosenPhenoParentData.galPhenotype === "Gal_ON" ? 0 : chosenPhenoParentData.gensSinceGalOn + 1,
            lagTimeRemaining: 0.0 // Reset lag, adapt below
        };

        const G_at_spot = glucose.get(spotCoord) || 0; const H_at_spot = galactose.get(spotCoord) || 0;
        let isGalPresentAtSpot = H_at_spot >= nutrientThresholdForPresence;
        let isGluHighAtSpot = G_at_spot >= glucoseThresholdForRepression;

        if (daughterCell.galPhenotype === "Gal_ON") {
            if (!isGalPresentAtSpot || isGluHighAtSpot) {
                daughterCell.galPhenotype = "Gal_OFF"; daughterCell.gensSinceGalOn = 1;
            } else {
                daughterCell.gensSinceGalOn = 0; // Still ON
            }
        } else { // galPhenotype === "Gal_OFF"
            if (isGalPresentAtSpot && !isGluHighAtSpot) {
                let memoryFactor = exp(-daughterCell.gensSinceGalOn / memoryDecayGens);
                let lagToEnter = minLagTimeWithMemory + (maxLagTime - minLagTimeWithMemory) * (1.0 - memoryFactor);
                if (lagTimeStochasticity > 0) lagToEnter += randomGaussian() * lagToEnter * lagTimeStochasticity;
                daughterCell.lagTimeRemaining = max(0.0, lagToEnter);
            }
        }
        nextLiveCells.set(spotCoord, daughterCell);
        newbornCoordsAndData.push({coordStr: spotCoord, cellData: daughterCell});
        successfulBirthsCount++;
        parentsUsedThisStep.add(phenoParentCoord);
        const [nq_b, nr_b] = spotCoord.split(',').map(Number); // Update colony radius with newborn
        currentColonyMaxRadius = max(currentColonyMaxRadius, axialDistance(nq_b, nr_b, 0, 0));
    }

    // --- PASS 5: Update Lag Times for ALL Cells (v30 logic) ---
    const elapsedTimeThisStep = successfulBirthsCount * dtPerCell;
    if (elapsedTimeThisStep > epsilon) {
        for (const [coordStr, cell] of nextLiveCells.entries()) {
            if (cell.lagTimeRemaining > epsilon) {
                cell.lagTimeRemaining -= elapsedTimeThisStep;
                if (cell.lagTimeRemaining <= epsilon) {
                    cell.lagTimeRemaining = 0.0; cell.galPhenotype = "Gal_ON"; cell.gensSinceGalOn = 0;
                }
            }
        }
    }

    // --- PASS 6: Update Frontier Set based on births (v30 logic) ---
    let finalNextFrontierSet = new Set();
    for (const {coordStr: newbornCoordStr} of newbornCoordsAndData) {
        finalNextFrontierSet.add(newbornCoordStr);
        const [q, r] = newbornCoordStr.split(',').map(Number);
        for (const dir of axialDirections) {
            const nq = q + dir.q; const nr = r + dir.r;
            const neighborCoordStr = `${nq},${nr}`;
            if (nextLiveCells.has(neighborCoordStr) && neighborCoordStr !== newbornCoordStr) {
                const [nq_n, nr_n] = neighborCoordStr.split(',').map(Number); let hasEmpty = false;
                for(const dir2 of axialDirections){ const nnq = nq_n + dir2.q; const nnr = nr_n + dir2.r; if(glucose.has(`${nnq},${nnr}`) && !nextLiveCells.has(`${nnq},${nnr}`)){ hasEmpty = true; break; } }
                if(hasEmpty) finalNextFrontierSet.add(neighborCoordStr);
            }
        }
    }
    // Refine frontier set
    let evenMoreFinalFrontierSet = new Set();
    for(const coord of finalNextFrontierSet) {
        if (nextLiveCells.has(coord)) {
             const [q_f, r_f] = coord.split(',').map(Number); let hasEmpty = false;
             for(const dir of axialDirections){ const nq_ef = q_f + dir.q; const nr_ef = r_f + dir.r; if(glucose.has(`${nq_ef},${nr_ef}`) && !nextLiveCells.has(`${nq_ef},${nr_ef}`)){ hasEmpty = true; break; } }
             if(hasEmpty) evenMoreFinalFrontierSet.add(coord);
        }
    }


    // --- PASS 7: Update Nutrient Profile if Dynamic Bands Mode ---
    if (PARAMS.gradientMode === "random_alternating_bands") {
        if (ensureNutrientBandsGenerated()) { // Returns true if bands were changed
            populateNutrientMapsAndDrawBuffer(); // Redraw buffer and update G/H maps
        }
    }

    // --- PASS 8: Final State Update ---
    liveCells = nextLiveCells;
    frontierSet = evenMoreFinalFrontierSet; // Use the most refined set
    lastStepBirths = successfulBirthsCount;
    if (liveCells.size === 0) lastStepBirths = 0;
}


// =============================================================================
// --- Drawing Functions (Cell drawing uses new colors) ---
// =============================================================================
function drawLiveCells(targetBuffer){
    if (!targetBuffer) return;
    for (const [coordStr, cell] of liveCells.entries()) {
         const [q, r] = coordStr.split(',').map(Number);
         let cellCol;
         if (cell.lagTimeRemaining > epsilon) {
             cellCol = cell.genotype === 1 ? cellColor_Green_Lag : cellColor_Red_Lag;
         } else if (cell.galPhenotype === "Gal_ON") {
             cellCol = cell.genotype === 1 ? cellColor_Green_GalON : cellColor_Red_GalON;
         } else { // Gal_OFF
             cellCol = cell.genotype === 1 ? cellColor_Green_GalOFF : cellColor_Red_GalOFF;
         }
         drawHexagon(targetBuffer, q, r, cellCol, color(50), 0.3); // Thinner, darker stroke
    }
}

function drawLiveCellsOnBuffer(targetBuffer){
    if (!targetBuffer) return;
    targetBuffer.push();
    targetBuffer.translate(targetBuffer.width / 2, targetBuffer.height / 2);
    targetBuffer.clear();
    drawLiveCells(targetBuffer);
    targetBuffer.pop();
 }

// =============================================================================
// --- UI Setup & Updates (SIMPLIFIED) ---
// =============================================================================
function createUI() {
    const uiPanel = select('#ui-panel'); uiPanel.html(''); // Clear previous UI
    uiElements = { controlGroups: {}, sliders: {}, selectors: {}, buttons: {} };
    uiParamDisplays = {};

    const createSliderRow = (label, paramKey, minVal, maxVal, stepVal, formatDigits = 2, subKey = null, isNutrientParam = false) => {
        let path = subKey ? PARAMS[paramKey][subKey] : PARAMS[paramKey];
        let rowDiv = createDiv(label).parent(uiPanel).addClass('ui-row');
        let slider = createSlider(minVal, maxVal, path, stepVal)
            .parent(rowDiv)
            .input(() => {
                let val = slider.value();
                if (subKey) PARAMS[paramKey][subKey] = val; else PARAMS[paramKey] = val;
                select('#' + (subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay').html(nf(val, 0, formatDigits));
                if (isNutrientParam && !start && PARAMS.gradientMode === "random_alternating_bands") {
                    console.log(`Nutrient profile param ${paramKey}${subKey ? '.'+subKey : ''} changed, re-init nutrients...`);
                    initGridAndNutrients(); redraw();
                } else if ((paramKey.includes('gradient') || paramKey.includes('noise') || paramKey.includes('stripe')) && !start && PARAMS.gradientMode !== "random_alternating_bands") {
                     initGridAndNutrients(); redraw();
                }
            });
        createSpan(nf(path, 0, formatDigits)).parent(rowDiv).id((subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay');
        uiElements.sliders[(subKey ? paramKey + '_' + subKey : paramKey)] = slider;
        uiParamDisplays[(subKey ? paramKey + '_' + subKey : paramKey)] = select('#' + (subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay');
        return rowDiv;
    };
    const createSectionHeader = (text) => createP(`<b>${text}</b>`).parent(uiPanel);
    const createSelectRow = (label, paramKey, options, isProfileRelated = false) => {
        let rowDiv = createDiv(label + ": ").parent(uiPanel).addClass('ui-row');
        let sel = createSelect().parent(rowDiv);
        options.forEach(opt => sel.option(typeof opt === 'string' ? opt : opt.label, typeof opt === 'string' ? opt : opt.value));
        sel.selected(PARAMS[paramKey]);
        sel.changed(() => {
            let prevMode = PARAMS[paramKey];
            PARAMS[paramKey] = sel.value();
            console.log(`${paramKey} changed to ${sel.value()}`);
            if ((paramKey === 'gradientMode' || isProfileRelated) && !start) {
                initGridAndNutrients(); redraw();
                if (paramKey === 'gradientMode') updateUIVisibility(); // Show/hide relevant old mode controls
            } else if (paramKey === 'initialConditionType' && !start) {
                 // No grid reset needed, just affects planting.
            }
        });
        uiElements.selectors[paramKey] = sel; return rowDiv;
    };

    createSectionHeader("Simulation Controls:");
    let controlsDiv = createDiv().parent(uiPanel).id('controls-container');
    uiElements.buttons.startBtn = createButton("Start / Pause (Space)").parent(controlsDiv).mousePressed(toggleSimulation);
    uiElements.buttons.resetBtn = createButton("Reset Sim").parent(controlsDiv).mousePressed(() => { initGridAndNutrients(); redraw(); });
    uiElements.buttons.gradientOverlayBtn = createButton("Show Gradient").parent(controlsDiv).mousePressed(toggleGradientOverlay);

    createSectionHeader("Core Parameters:");
    createSelectRow("Nutrient Mode:", 'gradientMode',
        ["random_alternating_bands", 'smooth_radial', 'smooth_linear_h', 'smooth_linear_v', 'gaussian_source', 'noisy', 'stripes_h', 'stripes_v', 'uniform_g', 'uniform_h', 'uniform_mix']
    );
    // Show these only if NOT random_alternating_bands
    uiElements.controlGroups.oldGradientControls = [];
    createSliderRow("Gradient Center (%):", 'gradientCenter', 0.0, 1.0, 0.01, 2).parent(uiElements.controlGroups.oldGradientControls[0] = createDiv().parent(uiPanel)); // Group old params
    createSliderRow("Transition Width (%):", 'gradientTransitionWidth', 0.05, 1.0, 0.01, 2).parent(uiElements.controlGroups.oldGradientControls[0]);

    createSelectRow("Initial Cell State:", 'initialConditionType', ["mix", "all_green", "all_red"]);
    createSliderRow("Initial Green Fraction (if mix):", 'initialGenotypeMix', 0, 1, 0.01, 2);
    createSliderRow("Max Lag Time:", 'maxLagTime', 0, 2.0, 0.01, 2);
    createSliderRow("Memory Decay (gens):", 'memoryDecayGens', 0.1, 10.0, 0.1, 1);

    createSectionHeader("Nutrient Profile (Random Bands):");
    const np = 'nutrientProfile'; // Alias
    createSelectRow("Homeland Env Type:", np, [{label: "Glu-High / Gal-Low", value:"GluHigh_GalLow"}, {label:"Glu-Low / Gal-High", value:"GluLow_GalHigh"}], true, 'homeland_env_type');
    createSliderRow("Homeland Radius (R0):", np, 0.1, 10.0, 0.1, 1, 'R0_homeland_radius', true);
    createSliderRow("Min Band Width:", np, 1.0, 10.0, 0.1, 1, 'min_band_width', true);
    createSliderRow("Max Band Width:", np, 1.0, 15.0, 0.1, 1, 'max_band_width', true);

    createSectionHeader("Status:");
    let statusDiv = createDiv().parent(uiPanel).id('status-section');
    createSpan('Step: 0').parent(statusDiv).id('stepValueDisplay');
    createSpan('Cell Count: 0').parent(statusDiv).id('cellCountValueDisplay');
    createSpan('Max Radius: 0').parent(statusDiv).id('maxRadiusValueDisplay');
    createSpan('Nutrient Bands: 0').parent(statusDiv).id('nutrientBandsCountDisplay');
    createSpan('Status: Initialized').parent(statusDiv).id('statusValueDisplay');

    updateUIVisibility(); // Initial call
    updateUIParamDisplays();
}

function updateUIVisibility() {
    const isRandomBandsMode = PARAMS.gradientMode === "random_alternating_bands";
    // Hide/show old gradient controls
    if (uiElements.controlGroups.oldGradientControls && uiElements.controlGroups.oldGradientControls[0]) {
        if (isRandomBandsMode) uiElements.controlGroups.oldGradientControls[0].hide();
        else uiElements.controlGroups.oldGradientControls[0].show();
    }
    // Could add more logic here for other old modes if their UI sections were re-added.
    // For now, "Nutrient Profile (Random Bands)" section is always visible,
    // but its parameters only have an effect if gradientMode is "random_alternating_bands".
}


function toggleSimulation() { start = !start; if (start) loop(); else noLoop(); updateDynamicUIText(); if (!start) redraw(); }
function toggleGradientOverlay() { showGradientOverlay = !showGradientOverlay; uiElements.buttons.gradientOverlayBtn.html(showGradientOverlay ? 'Show Cells' : 'Show Gradient'); if (!start) redraw(); }

function updateUIParamDisplays() { // Simplified: just update values, visibility handled by updateUIVisibility
    for (const key in uiParamDisplays) {
        let valueToDisplay, pKey = key, subKey = null, digits = 2;
        if (key.includes('_')) {
            const parts = key.split('_'); pKey = parts[0]; subKey = parts.slice(1).join('_');
            if (PARAMS[pKey] && typeof PARAMS[pKey] === 'object' && PARAMS[pKey].hasOwnProperty(subKey)) valueToDisplay = PARAMS[pKey][subKey]; else continue;
        } else if (PARAMS.hasOwnProperty(key)) valueToDisplay = PARAMS[key]; else continue;

        // Infer digits (very simplified for this version)
        if (key.includes('Rate')) digits = 5;
        else if (key.includes('Radius') || key.includes('Width') || key.includes('Decay')) digits = 1;
        else if (key.includes('initialPatchRadius')) digits = 0;

        if (uiParamDisplays[key] && uiParamDisplays[key].html() !== nf(valueToDisplay, 0, digits)) {
            uiParamDisplays[key].html(nf(valueToDisplay, 0, digits));
        }
    }
}

function updateDynamicUIText() {
    select('#stepValueDisplay')?.html(`Step: ${simulationSteps}`);
    select('#cellCountValueDisplay')?.html(`Cell Count: ${liveCells.size}`);
    select('#maxRadiusValueDisplay')?.html(`Max Radius: ${currentColonyMaxRadius.toFixed(1)}`);
    select('#nutrientBandsCountDisplay')?.html(`Nutrient Bands: ${currentNutrientBands.length + (PARAMS.nutrientProfile.R0_homeland_radius > 0 ? 1:0)}`);
    let statusSpan = select('#statusValueDisplay');
    if (statusSpan) {
        let statusMsg = "Initialized", statusColor = "#6c757d";
        if (start) { statusMsg = "RUNNING"; statusColor = "#28a745"; }
        else if (simulationSteps > 0) {
            if (liveCells.size === 0) { statusMsg = "EXTINCTION"; statusColor = "#dc3545"; }
            else if (lastStepBirths === 0 && finalStateRendered) { statusMsg = "STASIS - Stopped"; statusColor = "#dc3545"; }
            else if (lastStepBirths === 0) { statusMsg = "PAUSED (Stasis likely)"; statusColor = "#ffc107"; }
            else { statusMsg = "PAUSED"; statusColor = "#ffc107"; }
        } else { statusMsg = "Ready / PAUSED"; statusColor = "#6c757d"; }
        statusSpan.html(`Status: ${statusMsg}`).style('color', statusColor);
    }
}

// =============================================================================
// --- User Interaction (Planting uses new cell defaults) ---
// =============================================================================
function mousePressed() {
<<<<<<< HEAD
    if (mouseX >= 0 && mouseX <= width && mouseY >= 0 && mouseY <= height) {
        if (mouseX > canvasWidth) return; // Ignore clicks in UI panel
        let clickX_in_grid_space = mouseX - gridOriginX; let clickY_in_grid_space = mouseY - gridOriginY;
        const hexCoords = pixelToHex(clickX_in_grid_space, clickY_in_grid_space);
        if (axialDistance(hexCoords.q, hexCoords.r, 0, 0) <= maxWorldRadius + 2) {
            plantPatch(hexCoords.q, hexCoords.r, PARAMS.initialPatchRadius);
            if (!start) redraw();
        } else console.log("Click outside defined world radius.");
    }
}

function plantPatch(centerQ, centerR, radius) {
    let count = 0; let newlyPlantedCoords = [];
    const { initialConditionType, initialGenotypeMix, memoryDecayGens } = PARAMS;

    for (let q = -radius; q <= radius; q++) {
        let r1 = Math.max(-radius, -q - radius); let r2 = Math.min(radius, -q + radius);
        for (let r = r1; r <= r2; r++) {
            let targetQ = centerQ + q; let targetR = centerR + r;
            const coordStr = `${targetQ},${targetR}`;
            if (glucose.has(coordStr) && !liveCells.has(coordStr)) {
                let genotype_val = 0; // Default Red
                if (initialConditionType === "all_green") genotype_val = 1;
                else if (initialConditionType === "mix" && random() < initialGenotypeMix) genotype_val = 1;
                const newCell = { genotype: genotype_val, galPhenotype: "Gal_OFF", gensSinceGalOn: memoryDecayGens * 3, lagTimeRemaining: 0.0 };
                liveCells.set(coordStr, newCell); newlyPlantedCoords.push(coordStr); count++;
                currentColonyMaxRadius = max(currentColonyMaxRadius, axialDistance(targetQ, targetR, 0,0));
            }
        }
    }
=======
    let canvasRect = select('#canvas-container').elt.getBoundingClientRect();
    let bodyRect = document.body.getBoundingClientRect();
    let clickXInCanvas = mouseX - (canvasRect.left - bodyRect.left);
    let clickYInCanvas = mouseY - (canvasRect.top - bodyRect.top);

    if (clickXInCanvas >= 0 && clickXInCanvas <= canvasWidth && clickYInCanvas >= 0 && clickYInCanvas <= height) {
        let clickXRelativeToGridOrigin = clickXInCanvas - gridOriginX; let clickYRelativeToGridOrigin = clickYInCanvas - gridOriginY;
        const hexCoords = pixelToHex(clickXRelativeToGridOrigin, clickYRelativeToGridOrigin);
        if (axialDistance(hexCoords.q, hexCoords.r, 0, 0) <= maxInitRadius + 2) { plantPatch(hexCoords.q, hexCoords.r, 2); if (!start) redraw(); }
        else console.log("Click outside defined grid radius.");
    }
}

function plantPatch(centerQ, centerR, radius) {
    let pValueToPlant = constrain(PARAMS.initialP, 0, 1); let count = 0; let newlyPlantedCoords = [];
    for (let q = -radius; q <= radius; q++) { let r1 = Math.max(-radius, -q - radius); let r2 = Math.min(radius, -q + radius); for (let r = r1; r <= r2; r++) { let targetQ = centerQ + q; let targetR = centerR + r; const coordStr = `${targetQ},${targetR}`; if (glucose.has(coordStr) && !liveCells.has(coordStr)) { const newCell = { p: pValueToPlant }; liveCells.set(coordStr, newCell); newlyPlantedCoords.push(coordStr); count++; } } }
>>>>>>> parent of 774b6ad (Update sketch.js)
    if (count > 0) {
        let cellsToConsiderForFrontier = new Set(newlyPlantedCoords);
        for (const coordStr of newlyPlantedCoords) {
            const [q_p, r_p] = coordStr.split(',').map(Number);
            for (const dir of axialDirections) {
                const nq = q_p + dir.q; const nr = r_p + dir.r;
                const neighborCoordStr = `${nq},${nr}`;
                if (liveCells.has(neighborCoordStr) && !newlyPlantedCoords.includes(neighborCoordStr)) cellsToConsiderForFrontier.add(neighborCoordStr);
            }
        }
        for (const coordStrToEvaluate of cellsToConsiderForFrontier) {
            const [q_e, r_e] = coordStrToEvaluate.split(',').map(Number); let hasEmptyNeighbor = false;
            for (const dir of axialDirections) { const nq = q_e + dir.q; const nr = r_e + dir.r; if (glucose.has(`${nq},${nr}`) && !liveCells.has(`${nq},${nr}`)) { hasEmptyNeighbor = true; break; } }
            if (hasEmptyNeighbor) frontierSet.add(coordStrToEvaluate); else frontierSet.delete(coordStrToEvaluate);
        }
        console.log(`Planted ${count} cells @ (${centerQ},${centerR}). Init: ${initialConditionType}. Frontier updated.`);
        lastStepBirths = count; if (!start && finalStateRendered) finalStateRendered = false;
        if (!start) { updateDynamicUIText(); redraw(); }
        // If new cells expand colony radius significantly, may need to update nutrient bands now
        if (PARAMS.gradientMode === "random_alternating_bands") {
            if (ensureNutrientBandsGenerated()) populateNutrientMapsAndDrawBuffer();
        }
    } else console.log(`Could not plant @ (${centerQ},${centerR}) - occupied or outside range.`);
}
<<<<<<< HEAD

function keyPressed() { if (key === ' ') { toggleSimulation(); return false; } }
=======
>>>>>>> parent of 774b6ad (Update sketch.js)
// =============================================================================
