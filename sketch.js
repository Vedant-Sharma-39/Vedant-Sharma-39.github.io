// --- Core Microbial Strategy Simulation (v30 - Julia Model Port: Genotype, Phenotype, Lag, Dynamic Bands) ---
// Features:
// - Ported core logic from a Julia model:
//   - Discrete Genotypes (Red/Green) instead of continuous 'p'.
//   - Galactose Phenotype (Gal_ON / Gal_OFF) with switching.
//   - Lag phase for Gal_OFF -> Gal_ON switching, influenced by memory.
//   - Dynamic Banded Radial Nutrient Environment: Bands of high Glc/Gal can change.
//   - Specific fitness rules based on genotype, phenotype, and nutrients.
// - Colonization Model: Viability Threshold, Collective Pressure, Probabilistic Success, Parent Limitation.
// - User interaction: Click to plant, Spacebar Start/Pause, Gradient Overlay Toggle.
// --- NO CELL DEATH in this version (inherited from v29) ---
// ---

// --- Simulation Settings ---
let maxWorldRadius = 50; // Renamed from maxInitRadius to clarify it's the world boundary
let hexSize = 5;
const epsilon = 1e-9;
let drawInterval = 2; // Keep at 1 for smoother dynamic gradient updates if performance allows

// --- MODEL PARAMETERS (Defaults, adjustable via UI) ---
let PARAMS = {
  // --- Cell Genotype, Phenotype & Lag (NEW from Julia Model) ---
  initialGenotypeMix: 0.5, // Fraction of Green (1) cells if initial_condition_type is 'mix'
  initialConditionType: "mix", // "mix", "all_green", "all_red"
  mutationRateGreenToRed: 0.001, // mu_f_genotype (Green=1 -> Red=0)
  mutationRateRedToGreen: 0.0001, // mu_b_genotype (Red=0 -> Green=1)

  maxLagTime: 0.5,      // Max lag time (fraction of a "generation", time units arbitrary)
  memoryDecayGens: 3.0, // Generations for memory to decay significantly
  minLagTimeWithMemory: 0.05, // Min lag time if cell has recent memory of Gal
  lagTimeStochasticity: 0.1, // StdDev fraction for lag time noise (0 for no noise)
  dtPerCell: 0.025,     // Effective time step advanced per cell added (for lag reduction)

  // --- Fitness Parameters (NEW from Julia Model) ---
  fitnessParams: {
    // Fitness = Growth Rate relative to a baseline (e.g., 1.0 = 1 division/generation)
    // Green Genotype (1), Red Genotype (0)
    fitness_green_on_gal: 1.5,  // Green, Gal_ON, consuming Gal
    fitness_red_on_gal: 0.8,    // Red,   Gal_ON, consuming Gal
    fitness_on_no_gal_cost_factor: 0.7, // Multiplier for Glc fitness if Gal_ON but no/low Gal (or Glc repressed)
    fitness_green_off_glu: 1.0, // Green, Gal_OFF, consuming Glu
    fitness_red_off_glu: 0.95,  // Red,   Gal_OFF, consuming Glu
    fitness_off_on_gal: 0.05,   // Gal_OFF, but Gal present (and no Glc repression) - low fitness, should be inducing
    fitness_during_lag: 0.01,   // Fitness while in lag phase
  },
  // Nutrient thresholds for phenotype switching & fitness logic
  nutrientThresholdForPresence: 0.5, // Min nutrient level to be considered "present" for phenotype switching
  glucoseThresholdForRepression: 0.7, // High Glc level that represses Gal even if Gal present

  // --- Colonization (Largely from v29, K_G, K_H, h_max removed) ---
  minColonizationFitness: 0.05, // Min fitness for a cell to contribute to colonization pressure
  colonizationPressureFactor: 0.8,
  maxColonizationProbability: 0.95,
  initialPatchRadius: 2, // Radius of patch planted on click

  // --- Nutrient Gradient Core ---
  gradientMode: "dynamic_banded_radial", // Defaulting to new mode

  // --- Dynamic Banded Radial Controls (NEW) ---
  nutrientBandDefaultGlu: 0.1,
  nutrientBandDefaultGal: 0.1,
  nutrientBandHighGlu: 1.0,
  nutrientBandHighGal: 1.0,
  nutrientBandChangeIntervalCells: 150, // Update bands every X cells added
  nutrientBandMaxRadiusFactor: 1.2, // Bands can extend up to current colony radius * this factor
  nutrientBandWidthMin: 3.0,
  nutrientBandWidthMax: 8.0,
  nutrientBandMinStartRadius: 1.0, // Min distance from origin for a band to start
  nutrientBandMinSeparation: 2.0, // Min separation between centers of different nutrient bands

  // --- Smooth Radial/Linear Controls (from v29, for other modes) ---
  gradientCenter: 0.25,
  gradientTransitionWidth: 0.5,

  // --- Gaussian Source Controls (from v29) ---
  gaussianSourceQ: 0, gaussianSourceR: 0, gaussianWidth: 15,
  gaussianPeakG: 1.0, gaussianPeakH: 0.0, gaussianBaseG: 0.0, gaussianBaseH: 0.2,

  // --- Noisy Gradient Controls (from v29) ---
  noiseScale: 0.1, noiseIntensityG: 0.3, noiseIntensityH: 0.3,
  noiseBaseGradient: "uniform_mix",

  // --- Stripe Controls (from v29) ---
  stripeWidth: 40,
};

// --- Simulation State ---
let glucose = new Map(); let galactose = new Map(); let liveCells = new Map();
let frontierSet = new Set(); let start = false; let simulationSteps = 0; let lastStepBirths = 1;
let showGradientOverlay = false;
let currentColonyMaxRadius = 0; // For dynamic nutrient bands
let cellsAddedSinceNutrientChange = 0; // For dynamic nutrient bands
// Store actual band locations: array of {r_in: val, r_out: val} objects
let currentNutrientGluBands = [];
let currentNutrientGalBands = [];


// --- Graphics & Grid ---
let nutrientBuffer, finalStateBuffer; let finalStateRendered = false; const SQRT3 = Math.sqrt(3); let hexWidth, hexHeight, vertSpacing;
const axialDirections = [ { q: +1, r: 0 }, { q: +1, r: -1 }, { q: 0, r: -1 }, { q: -1, r: 0 }, { q: -1, r: +1 }, { q: 0, r: +1 } ];
let canvasWidth, canvasHeight, bufferWidth, bufferHeight, gridOriginX, gridOriginY;

// --- Color Definitions ---
let nutrientColorG, nutrientColorH, nutrientColorMix;
// New cell colors based on Genotype (0=Red, 1=Green) and Phenotype (Gal_ON/Gal_OFF/Lag)
let cellColor_Green_GalOFF, cellColor_Green_GalON, cellColor_Green_Lag;
let cellColor_Red_GalOFF, cellColor_Red_GalON, cellColor_Red_Lag;


// --- UI Elements ---
let uiElements = { controlGroups: {} };
let uiParamDisplays = {};

// =============================================================================
// --- Math & Grid Helpers (Mostly unchanged from v29) ---
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
// --- Drawing Helpers (Mostly unchanged) ---
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
// =================================================S============================
function setup() {
    calculateGridDimensions();
    nutrientBuffer = createGraphics(bufferWidth, bufferHeight); finalStateBuffer = createGraphics(bufferWidth, bufferHeight); finalStateBuffer.pixelDensity(1);
    console.log(`Grid/Buffer Size: ${bufferWidth}x${bufferHeight}. Hex Size: ${hexSize}`);

    // Nutrient Colors (for gradient display)
    nutrientColorG = color(100, 220, 100, 200); nutrientColorH = color(200, 100, 220, 200); nutrientColorMix = color(240);
    // Cell Colors (Genotype & Phenotype)
    cellColor_Green_GalOFF = color(144, 238, 144); // lightgreen
    cellColor_Green_GalON  = color(0, 100, 0);     // darkgreen
    cellColor_Green_Lag    = color(255, 255, 0);   // yellow
    cellColor_Red_GalOFF   = color(250, 128, 114); // salmon
    cellColor_Red_GalON    = color(139, 0, 0);      // darkred
    cellColor_Red_Lag      = color(255, 165, 0);   // orange


    let totalCanvasHeight = Math.max(canvasHeight, 1100); // Ensure height for UI panel
    let mainCanvas = createCanvas(canvasWidth, totalCanvasHeight);
    let canvasContainer = select('#canvas-container'); if (canvasContainer) mainCanvas.parent('canvas-container'); else console.error("#canvas-container not found.");
    noiseDetail(4, 0.5);
    createUI(); // UI creation needs to happen after PARAMS are defined
    initGridAndNutrients(); // Combined initialization
    background(240); push(); translate(gridOriginX, gridOriginY); image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2); drawLiveCells(this); pop();
    updateUIParamDisplays(); updateDynamicUIText(); noLoop();
}

function draw() {
     if (start) { stepSimulation(); }
    // Redraw overlay or cells more frequently if needed, esp. if nutrients change dynamically
    if (frameCount % drawInterval === 0 || !start || showGradientOverlay || (PARAMS.gradientMode === "dynamic_banded_radial" && start)) {
        clear(); push(); translate(gridOriginX, gridOriginY);
        image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2); // Nutrient buffer might have been updated
        if (!showGradientOverlay) {
             if (start) { drawLiveCells(this); }
             else {
                 if (finalStateRendered) image(finalStateBuffer, -bufferWidth / 2, -bufferHeight / 2);
                 else drawLiveCells(this);
             }
        }
        pop();
    }
    updateDynamicUIText(); updateUIParamDisplays(); // Keep UI responsive
}

// =============================================================================
// --- Nutrient Environment Logic (NEW / HEAVILY MODIFIED) ---
// =============================================================================

// Gets G/H levels for a specific hex (q,r) based on the CURRENT nutrient settings
function getNutrientLevelsAtHex(q, r) {
    const { gradientMode, noiseScale, noiseIntensityG, noiseIntensityH, noiseBaseGradient,
            gaussianSourceQ, gaussianSourceR, gaussianWidth, gaussianPeakG, gaussianPeakH, gaussianBaseG, gaussianBaseH,
            stripeWidth, gradientCenter, gradientTransitionWidth,
            nutrientBandDefaultGlu, nutrientBandDefaultGal, nutrientBandHighGlu, nutrientBandHighGal
          } = PARAMS;

    let gVal = 0.5, hVal = 0.5; // Default for uniform_mix or fallback

    const R_world = maxWorldRadius; // Max extent of the grid
    const pixelPos = hexToPixel(q, r); // For linear/stripe gradients

    // Calculate effective min/max pixel coordinates for normalization if needed by a mode
    // This is a bit of a simplification; ideally, pre-calculate once for the whole grid if static.
    // For dynamic bands, this is not directly used for G/H values but for other modes.
    const approxMinPixelX = -R_world * hexWidth / 2;
    const approxMaxPixelX = R_world * hexWidth / 2;
    const approxMinPixelY = -R_world * vertSpacing / 2;
    const approxMaxPixelY = R_world * vertSpacing / 2;

    if (gradientMode === "dynamic_banded_radial") {
        let dist = axialDistance(q, r, 0, 0);
        gVal = nutrientBandDefaultGlu;
        hVal = nutrientBandDefaultGal;
        for (const band of currentNutrientGluBands) {
            if (dist > band.r_in && dist <= band.r_out) {
                gVal = nutrientBandHighGlu;
                break;
            }
        }
        for (const band of currentNutrientGalBands) {
            if (dist > band.r_in && dist <= band.r_out) {
                hVal = nutrientBandHighGal;
                break;
            }
        }
    } else { // Fallback to v29 static gradient logic
        let baseValues = calculateBaseGradientValue_Static( // Use the old helper for static modes
            gradientMode === 'noisy' ? noiseBaseGradient : gradientMode,
            q, r, R_world, pixelPos, approxMinPixelX, approxMaxPixelX, approxMinPixelY, approxMaxPixelY
        );
        gVal = baseValues.g; hVal = baseValues.h;

        if (gradientMode === 'noisy') {
            let noiseValG = noise(q * noiseScale, r * noiseScale, 0.1);
            let noiseValH = noise(q * noiseScale + 100, r * noiseScale + 100, 0.2);
            gVal += (noiseValG - 0.5) * 2 * noiseIntensityG;
            hVal += (noiseValH - 0.5) * 2 * noiseIntensityH;
        }
    }
    return { g: constrain(gVal, 0, 1), h: constrain(hVal, 0, 1) };
}

// Renamed from calculateBaseGradientValue to avoid conflict, this is for STATIC modes from v29
function calculateBaseGradientValue_Static(mode, q, r, R, pixelPos, minPixelX, maxPixelX, minPixelY, maxPixelY) {
    const { gradientCenter, gradientTransitionWidth,
            gaussianSourceQ, gaussianSourceR, gaussianWidth, gaussianPeakG, gaussianPeakH, gaussianBaseG, gaussianBaseH,
            stripeWidth
          } = PARAMS;
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
    return { g: constrain(baseG, 0, 1), h: constrain(baseH, 0, 1) };
}


function initializeNutrientBands() {
    const { nutrientBandMinStartRadius, nutrientBandWidthMin, nutrientBandWidthMax } = PARAMS;
    // For initial setup, let's try one Glu band and one Gal band if space allows
    // Or, make it simpler: start with no specific bands, just defaults
    currentNutrientGluBands = [];
    currentNutrientGalBands = [];
    cellsAddedSinceNutrientChange = 0; // Reset counter

    // Example: Add one initial glucose band
    // let r_start_glu = random(nutrientBandMinStartRadius, maxWorldRadius / 3);
    // let w_glu = random(nutrientBandWidthMin, nutrientBandWidthMax);
    // if (r_start_glu + w_glu < maxWorldRadius * 0.8) { // Ensure it's within reasonable bounds
    //    currentNutrientGluBands.push({ r_in: r_start_glu, r_out: r_start_glu + w_glu });
    // }
    // console.log("Initial nutrient bands set (example: empty or one band).");
    // For now, let's call the dynamic update function to create initial bands
    updateNutrientBands(true); // true for initial setup
}

function updateNutrientBands(isInitialSetup = false) {
    const { nutrientBandMinStartRadius, nutrientBandWidthMin, nutrientBandWidthMax,
            nutrientBandMaxRadiusFactor, nutrientBandMinSeparation } = PARAMS;

    let effectiveMaxRadiusForBands = max(maxWorldRadius * 0.8, currentColonyMaxRadius * nutrientBandMaxRadiusFactor);
    effectiveMaxRadiusForBands = max(effectiveMaxRadiusForBands, nutrientBandMinStartRadius + nutrientBandWidthMax + 5); // ensure some space

    currentNutrientGluBands = [];
    currentNutrientGalBands = [];

    // Create one Glucose band
    if (effectiveMaxRadiusForBands > nutrientBandMinStartRadius + nutrientBandWidthMax) {
        let r_start_glu = random(nutrientBandMinStartRadius, effectiveMaxRadiusForBands - nutrientBandWidthMax);
        let w_glu = random(nutrientBandWidthMin, nutrientBandWidthMax);
        currentNutrientGluBands.push({ r_in: r_start_glu, r_out: r_start_glu + w_glu });
    }

    // Create one Galactose band, ensuring some separation from the Glucose band
    if (effectiveMaxRadiusForBands > nutrientBandMinStartRadius + nutrientBandWidthMax && currentNutrientGluBands.length > 0) {
        let attempts = 0;
        while (attempts < 20) {
            let r_start_gal = random(nutrientBandMinStartRadius, effectiveMaxRadiusForBands - nutrientBandWidthMax);
            let w_gal = random(nutrientBandWidthMin, nutrientBandWidthMax);

            let glu_center = currentNutrientGluBands[0].r_in + (currentNutrientGluBands[0].r_out - currentNutrientGluBands[0].r_in) / 2;
            let gal_center = r_start_gal + w_gal / 2;
            let dist_centers = abs(gal_center - glu_center);
            let min_dist_needed = (currentNutrientGluBands[0].r_out - currentNutrientGluBands[0].r_in) / 2 + w_gal / 2 + nutrientBandMinSeparation;

            if (dist_centers >= min_dist_needed) {
                currentNutrientGalBands.push({ r_in: r_start_gal, r_out: r_start_gal + w_gal });
                break;
            }
            attempts++;
        }
        if (attempts >= 20 && currentNutrientGalBands.length === 0 && effectiveMaxRadiusForBands > nutrientBandMinStartRadius + nutrientBandWidthMax * 2 + nutrientBandMinSeparation) {
            // Fallback if separation fails: place it far from the first band if possible
            let r_start_gal = (currentNutrientGluBands[0].r_out + nutrientBandMinSeparation + nutrientBandWidthMin < effectiveMaxRadiusForBands - nutrientBandWidthMin) ?
                              currentNutrientGluBands[0].r_out + nutrientBandMinSeparation : nutrientBandMinStartRadius;
            if (r_start_gal + nutrientBandWidthMin < effectiveMaxRadiusForBands) {
                 currentNutrientGalBands.push({ r_in: r_start_gal, r_out: r_start_gal + random(nutrientBandWidthMin, nutrientBandWidthMax) });
            }
        }
    } else if (effectiveMaxRadiusForBands > nutrientBandMinStartRadius + nutrientBandWidthMax && currentNutrientGluBands.length === 0) {
        // If no glu band, just place a gal band
        let r_start_gal = random(nutrientBandMinStartRadius, effectiveMaxRadiusForBands - nutrientBandWidthMax);
        let w_gal = random(nutrientBandWidthMin, nutrientBandWidthMax);
        currentNutrientGalBands.push({ r_in: r_start_gal, r_out: r_start_gal + w_gal });
    }


    cellsAddedSinceNutrientChange = 0;
    if (!isInitialSetup) {
        console.log(`Nutrient bands updated. MaxColR: ${currentColonyMaxRadius.toFixed(1)}. MaxBandR: ${effectiveMaxRadiusForBands.toFixed(1)} Glu: ${JSON.stringify(currentNutrientGluBands)}, Gal: ${JSON.stringify(currentNutrientGalBands)}`);
        // Redraw nutrient buffer as bands have changed
        populateNutrientMapsAndDrawBuffer(); // This will redraw the buffer
    }
}

function populateNutrientMapsAndDrawBuffer() {
    const R = maxWorldRadius;
    const bufferCenterX = bufferWidth / 2;
    const bufferCenterY = bufferHeight / 2;

    nutrientBuffer.push();
    nutrientBuffer.translate(bufferCenterX, bufferCenterY);
    nutrientBuffer.background(nutrientColorMix); // Base color for areas with no specific G/H dominance

    let bufferNutrientG_col = color(red(nutrientColorG), green(nutrientColorG), blue(nutrientColorG));
    let bufferNutrientH_col = color(red(nutrientColorH), green(nutrientColorH), blue(nutrientColorH));

    for (let q_ax = -R; q_ax <= R; q_ax++) {
        let r1 = Math.max(-R, -q_ax - R);
        let r2 = Math.min(R, -q_ax + R);
        for (let r_ax = r1; r_ax <= r2; r_ax++) {
            const coordStr = `${q_ax},${r_ax}`;
            const { g, h } = getNutrientLevelsAtHex(q_ax, r_ax);
            glucose.set(coordStr, g);
            galactose.set(coordStr, h);

            let totalNutrient = g + h;
            let hexColor;
            if (totalNutrient < epsilon) { // Very low nutrients
                hexColor = nutrientColorMix;
            } else if (g > h && g > epsilon) { // More G than H
                hexColor = lerpColor(nutrientColorMix, bufferNutrientG_col, constrain(g / (g + h + epsilon) * (g / 0.7),0,1) ); // Emphasize G if dominant
            } else if (h >= g && h > epsilon) { // More H than G or equal
                hexColor = lerpColor(nutrientColorMix, bufferNutrientH_col, constrain(h / (g + h + epsilon) * (h / 0.7),0,1) ); // Emphasize H if dominant
            } else { // Both very low or zero, but totalNutrient not < epsilon (should be covered by first case)
                hexColor = nutrientColorMix;
            }
            drawHexagon(nutrientBuffer, q_ax, r_ax, hexColor, null);
        }
    }
    nutrientBuffer.pop();
    if (start) redraw(); // Force redraw of main canvas to show updated buffer
}


// =============================================================================
// --- Initialization & Reset ---
// =============================================================================
function initGridAndNutrients() { // Combined and renamed
    start = false; loop(); lastStepBirths = 1; finalStateRendered = false; showGradientOverlay = false;
    if(uiElements.gradientOverlayBtn) uiElements.gradientOverlayBtn.html('Show Gradient');
    if(finalStateBuffer) finalStateBuffer.clear(); simulationSteps = 0;
    glucose.clear(); galactose.clear(); liveCells.clear(); frontierSet.clear();
    currentColonyMaxRadius = 0;
    cellsAddedSinceNutrientChange = 0;
    currentNutrientGluBands = [];
    currentNutrientGalBands = [];


    if (PARAMS.gradientMode === "dynamic_banded_radial") {
        initializeNutrientBands(); // Sets up currentNutrientGluBands, currentNutrientGalBands
    }
    // This populates glucose/galactose maps and draws the nutrientBuffer
    populateNutrientMapsAndDrawBuffer();

    console.log(`Grid Initialized. Grad: ${PARAMS.gradientMode}. Frontier Cleared.`);
    noLoop();
    updateDynamicUIText(); // Update status text
}

// =============================================================================
// --- Cell Fitness Logic (NEW from Julia Model) ---
// =============================================================================
function calculateCellFitness(cellData, G, H, fitnessP) {
    if (cellData.lagTimeRemaining > epsilon) {
        return fitnessP.fitness_during_lag;
    }

    const { nutrientThresholdForPresence, glucoseThresholdForRepression } = PARAMS;
    let isGalPresent = H > nutrientThresholdForPresence;
    // For catabolite repression, glucose needs to be significantly high
    let isGluHighForRepression = G > glucoseThresholdForRepression;

    if (cellData.galPhenotype === "Gal_ON") {
        if (isGalPresent && !isGluHighForRepression) { // Preferentially use Gal if ON and no strong Glu repression
            return cellData.genotype === 1 ? fitnessP.fitness_green_on_gal : fitnessP.fitness_red_on_gal;
        } else { // Gal-ON, but either no Gal or Glu is present & repressive (costly or using Glu)
            // Assume it uses Glucose if available, otherwise fitness is low
            let base_glu_fitness = G > nutrientThresholdForPresence ?
                                   (cellData.genotype === 1 ? fitnessP.fitness_green_off_glu : fitnessP.fitness_red_off_glu)
                                   : 0.01; // Very low fitness if Gal_ON and no utilizable nutrient
            return base_glu_fitness * fitnessP.fitness_on_no_gal_cost_factor;
        }
    } else { // Gal_OFF
        if (isGalPresent && !isGluHighForRepression) { // Should be inducing, but currently OFF
            return fitnessP.fitness_off_on_gal;
        } else { // Gal-OFF, and either no Gal or Glu is present (normal Glu metabolism or low nutrient)
            // Assume it uses Glucose if available
             return G > nutrientThresholdForPresence ?
                    (cellData.genotype === 1 ? fitnessP.fitness_green_off_glu : fitnessP.fitness_red_off_glu)
                    : 0.01; // Very low fitness if Gal_OFF and no Glucose
        }
    }
}


// =============================================================================
// --- Simulation Step Logic (HEAVILY MODIFIED) ---
// =============================================================================
function stepSimulation() {
    if (liveCells.size === 0 || (lastStepBirths === 0 && simulationSteps > 0)) {
        if (start) {
            let reason = liveCells.size === 0 ? "Extinction" : "Stasis";
            console.log(`${reason} @ step ${simulationSteps}. Stopping.`);
            if (liveCells.size > 0 && !finalStateRendered) {
                console.log("Rendering final state...");
                drawLiveCellsOnBuffer(finalStateBuffer);
                finalStateRendered = true;
                console.log("Final state rendered.");
            }
            start = false;
            noLoop();
            updateDynamicUIText();
        }
        return;
    }
    simulationSteps++;

    const { minColonizationFitness, colonizationPressureFactor, maxColonizationProbability,
            mutationRateGreenToRed, mutationRateRedToGreen,
            maxLagTime, memoryDecayGens, minLagTimeWithMemory, lagTimeStochasticity,
            dtPerCell, fitnessParams, nutrientThresholdForPresence, glucoseThresholdForRepression,
            nutrientBandChangeIntervalCells
          } = PARAMS;

    let nextLiveCells = new Map(liveCells);
    let nextFrontierSet = new Set();
    let newbornCoordsAndData = []; // Store {coordStr, cellData} for newborns
    let successfulBirthsCount = 0;

    // --- PASS 1: Calculate Fitness for all Current Frontier Cells ---
    let frontierFitness = new Map();
    currentColonyMaxRadius = 0; // Recalculate each step based on live cells
    for (const [coordStr, cell] of liveCells) { // Iterate all live cells to find max radius
        const [q,r] = coordStr.split(',').map(Number);
        currentColonyMaxRadius = max(currentColonyMaxRadius, axialDistance(q,r,0,0));
    }

    for (const coordStr of frontierSet) {
        if (!liveCells.has(coordStr)) continue;
        const cell = liveCells.get(coordStr);
        const G = glucose.get(coordStr) || 0;
        const H = galactose.get(coordStr) || 0;
        let fitness = calculateCellFitness(cell, G, H, fitnessParams);
        frontierFitness.set(coordStr, fitness);

        // Check if it should remain a frontier cell
        const [q, r] = coordStr.split(',').map(Number);
        let hasEmpty = false;
        for(const dir of axialDirections){
            const nq = q + dir.q; const nr = r + dir.r;
            const nCoordStr = `${nq},${nr}`;
            if(glucose.has(nCoordStr) && !liveCells.has(nCoordStr)){ // Check against current live cells
                hasEmpty = true; break;
            }
        }
        if(hasEmpty) nextFrontierSet.add(coordStr);
    }

    // --- PASS 2: Calculate Colonization Pressure & Chance for Each Relevant Empty Spot ---
    let potentialColonizations = new Map();
    let emptySpotsChecked = new Set();
    for (const frontierCoordStr of frontierSet) {
        if (!frontierFitness.has(frontierCoordStr)) continue;
        const [q, r] = frontierCoordStr.split(',').map(Number);
        for (const dir of axialDirections) {
            const nq = q + dir.q; const nr = r + dir.r;
            const emptySpotCoord = `${nq},${nr}`;
            if (glucose.has(emptySpotCoord) && !liveCells.has(emptySpotCoord) && !emptySpotsChecked.has(emptySpotCoord)) {
                emptySpotsChecked.add(emptySpotCoord);
                let viableNeighbors = []; // Store { coord, fitness, cellData }
                let totalPressure = 0;
                const [esq, esr] = emptySpotCoord.split(',').map(Number);

                for (const dir2 of axialDirections) {
                    const nnq = esq + dir2.q; const nnr = esr + dir2.r;
                    const neighborCoord = `${nnq},${nnr}`;
                    if (liveCells.has(neighborCoord) && frontierFitness.has(neighborCoord)) {
                        const neighborFitness = frontierFitness.get(neighborCoord);
                        if (neighborFitness >= minColonizationFitness) {
                            let pressure = map(neighborFitness, minColonizationFitness, 1.0, 0.01, 1.0, true);
                            totalPressure += pressure;
                            viableNeighbors.push({
                                coord: neighborCoord,
                                fitness: neighborFitness,
                                cellData: liveCells.get(neighborCoord) // Get full parent data
                            });
                        }
                    }
                }
                if (viableNeighbors.length > 0) {
                    let P_spot_colonized = constrain(totalPressure * colonizationPressureFactor, 0, maxColonizationProbability);
                    potentialColonizations.set(emptySpotCoord, { P_colonize: P_spot_colonized, viableNeighbors: viableNeighbors });
                }
            }
        }
    }

    // --- PASS 3: Determine Actual Colonizations (Stochastic Check) ---
    let colonizationEvents = []; // Stores { spotCoord, viableNeighbors }
    for (const [spotCoord, data] of potentialColonizations.entries()) {
        if (random() < data.P_colonize) {
            colonizationEvents.push({ spotCoord: spotCoord, viableNeighbors: data.viableNeighbors });
        }
    }

    // --- PASS 4: Resolve Parent Conflicts, Determine Daughter Genotype/Phenotype, Confirm Births ---
    shuffle(colonizationEvents, true); // Randomize order to resolve parent conflicts fairly
    let parentsUsedThisStep = new Set();

    for (const event of colonizationEvents) {
        const spotCoord = event.spotCoord;
        const viableParents = event.viableNeighbors; // viableNeighbors has {coord, fitness, cellData}

        // --- Select ONE actual parent for phenotypic inheritance (and check if it can reproduce) ---
        let totalParentFitness = 0;
        viableParents.forEach(p => totalParentFitness += p.fitness);

        if (totalParentFitness < epsilon) continue; // No viable parent can reproduce

        let chosenPhenoParentData = null;
        let phenoParentCoord = null;

        // Weighted selection for phenotypic parent, checking if parent already reproduced
        let cumulativeFitness = 0;
        let randomThreshold = random(totalParentFitness);
        let foundPhenoParent = false;

        // Shuffle parents to give fair chance if multiple have high fitness
        let shuffledViableParents = shuffle([...viableParents]);

        for (const parent of shuffledViableParents) {
            if (!parentsUsedThisStep.has(parent.coord)) {
                 cumulativeFitness += parent.fitness;
                 if (randomThreshold <= cumulativeFitness) {
                    chosenPhenoParentData = parent.cellData;
                    phenoParentCoord = parent.coord;
                    foundPhenoParent = true;
                    break;
                 }
            }
        }
        // Fallback if initial randomThreshold misses due to used parents
        if (!foundPhenoParent) {
            for (const parent of shuffledViableParents) {
                if (!parentsUsedThisStep.has(parent.coord)) {
                    chosenPhenoParentData = parent.cellData;
                    phenoParentCoord = parent.coord;
                    foundPhenoParent = true;
                    break;
                }
            }
        }

        if (!foundPhenoParent) continue; // All potential parents already reproduced

        // --- Daughter Genotype (based on ALL viable parents' fitness contributions) ---
        let weightedSumGreenContrib = 0.0;
        let weightedSumRedContrib = 0.0;
        for (const p of viableParents) { // Use original viableParents for genotype calc
            if (p.cellData.genotype === 1) weightedSumGreenContrib += p.fitness;
            else weightedSumRedContrib += p.fitness;
        }
        let pG_genotype = (weightedSumGreenContrib + weightedSumRedContrib > epsilon) ?
                          (weightedSumGreenContrib / (weightedSumGreenContrib + weightedSumRedContrib)) :
                          (weightedSumGreenContrib > 0 ? 1.0 : 0.0);

        let daughterGenotype = random() < pG_genotype ? 1 : 0; // 1 for Green, 0 for Red
        if (daughterGenotype === 1 && random() < mutationRateGreenToRed) daughterGenotype = 0;
        if (daughterGenotype === 0 && random() < mutationRateRedToGreen) daughterGenotype = 1;

        // --- Daughter's Initial Phenotypic State (Inherited) ---
        let daughterCell = {
            genotype: daughterGenotype,
            galPhenotype: chosenPhenoParentData.galPhenotype,
            gensSinceGalOn: chosenPhenoParentData.gensSinceGalOn + 1, // Tentatively increment
            lagTimeRemaining: chosenPhenoParentData.lagTimeRemaining // Inherit lag state
        };

        // --- Daughter Cell Adapts Phenotype to its New Environment ---
        const G_at_spot = glucose.get(spotCoord) || 0;
        const H_at_spot = galactose.get(spotCoord) || 0;
        let isGalPresentAtSpot = H_at_spot > nutrientThresholdForPresence;
        let isGluHighAtSpot = G_at_spot > glucoseThresholdForRepression;

        if (daughterCell.lagTimeRemaining <= epsilon) { // Only adapt if not currently mid-lag
            if (daughterCell.galPhenotype === "Gal_ON") {
                if (!isGalPresentAtSpot || isGluHighAtSpot) { // Condition to switch OFF
                    daughterCell.galPhenotype = "Gal_OFF";
                    daughterCell.gensSinceGalOn = 1; // Reset memory counter
                } else {
                    daughterCell.gensSinceGalOn = 0; // Still ON and utilizing Gal
                }
            } else { // galPhenotype === "Gal_OFF"
                if (isGalPresentAtSpot && !isGluHighAtSpot) { // Condition to try to switch ON
                    let memoryFactor = exp(-daughterCell.gensSinceGalOn / memoryDecayGens);
                    let lagToEnter = minLagTimeWithMemory +
                                   (maxLagTime - minLagTimeWithMemory) * (1.0 - memoryFactor);
                    if (lagTimeStochasticity > 0) {
                        lagToEnter += randomGaussian() * lagToEnter * lagTimeStochasticity; // randomGaussian is 0-mean, 1-stddev
                    }
                    daughterCell.lagTimeRemaining = max(0.0, lagToEnter);
                    // Phenotype switches to ON *after* lag completes (handled in Pass 5)
                }
                // else: remains Gal_OFF, gensSinceGalOn continues to count up (implicitly by inheritance)
            }
        }

        // --- Confirm Birth ---
        nextLiveCells.set(spotCoord, daughterCell);
        newbornCoordsAndData.push({coordStr: spotCoord, cellData: daughterCell}); // Keep track of new cell data
        successfulBirthsCount++;
        parentsUsedThisStep.add(phenoParentCoord); // Mark parent as used

        // Update colony max radius with the newborn
        const [nq, nr] = spotCoord.split(',').map(Number);
        currentColonyMaxRadius = max(currentColonyMaxRadius, axialDistance(nq, nr, 0, 0));
    }

    // --- PASS 5: Update Lag Times for ALL Cells (including newborns) ---
    const elapsedTimeThisStep = successfulBirthsCount * dtPerCell;
    if (elapsedTimeThisStep > epsilon) {
        for (const [coordStr, cell] of nextLiveCells.entries()) {
            if (cell.lagTimeRemaining > epsilon) {
                cell.lagTimeRemaining -= elapsedTimeThisStep;
                if (cell.lagTimeRemaining <= epsilon) {
                    cell.lagTimeRemaining = 0.0;
                    cell.galPhenotype = "Gal_ON"; // Switch phenotype after lag
                    cell.gensSinceGalOn = 0;
                }
            } else if (cell.galPhenotype === "Gal_OFF") {
                 // If not in lag and OFF, increment generations since Gal_ON,
                 // but only if it wasn't just born (newborns inherit and then might enter lag).
                 // This is tricky because newborns might have just entered lag.
                 // Let's assume gensSinceGalOn is primarily managed by inheritance and reset on ON.
                 // The +1 during inheritance handles the "time passing" aspect for memory.
            }
        }
    }

    // --- PASS 6: Update Frontier Set based on births ---
    // Add newborns to frontier consideration, and their existing live neighbors
    for (const {coordStr: newbornCoordStr} of newbornCoordsAndData) {
        nextFrontierSet.add(newbornCoordStr); // Newborn is on the frontier
        const [q, r] = newbornCoordStr.split(',').map(Number);
        for (const dir of axialDirections) {
            const nq = q + dir.q; const nr = r + dir.r;
            const neighborCoordStr = `${nq},${nr}`;
            if (nextLiveCells.has(neighborCoordStr) && neighborCoordStr !== newbornCoordStr) {
                // Check if this existing neighbor now needs to be on the frontier
                 let hasEmpty = false;
                 const [nq_n, nr_n] = neighborCoordStr.split(',').map(Number);
                 for(const dir2 of axialDirections){
                    const nnq = nq_n + dir2.q; const nnr = nr_n + dir2.r;
                    if(glucose.has(`${nnq},${nnr}`) && !nextLiveCells.has(`${nnq},${nnr}`)){
                        hasEmpty = true; break;
                    }
                 }
                 if(hasEmpty) nextFrontierSet.add(neighborCoordStr);
                 else nextFrontierSet.delete(neighborCoordStr); // If it was frontier but now surrounded
            }
        }
    }
    // Final filter for frontier: ensure they are still in nextLiveCells and have empty neighbors
    let finalNextFrontierSet = new Set();
    for (const coordStr of nextFrontierSet) {
        if (nextLiveCells.has(coordStr)) {
            const [q,r] = coordStr.split(',').map(Number);
            let hasEmpty = false;
            for(const dir of axialDirections){
                const nq_f = q + dir.q; const nr_f = r + dir.r;
                const nCoordStr_f = `${nq_f},${nr_f}`;
                if(glucose.has(nCoordStr_f) && !nextLiveCells.has(nCoordStr_f)){
                    hasEmpty = true; break;
                }
            }
            if(hasEmpty) finalNextFrontierSet.add(coordStr);
        }
    }


    // --- PASS 7: Update Nutrient Profile if Dynamic Bands Mode ---
    if (PARAMS.gradientMode === "dynamic_banded_radial") {
        cellsAddedSinceNutrientChange += successfulBirthsCount;
        if (cellsAddedSinceNutrientChange >= nutrientBandChangeIntervalCells) {
            updateNutrientBands(); // This will also redraw nutrientBuffer
        }
    }

    // --- PASS 8: Final State Update ---
    liveCells = nextLiveCells;
    frontierSet = finalNextFrontierSet;
    lastStepBirths = successfulBirthsCount;
    if (liveCells.size === 0) lastStepBirths = 0; // Ensure stasis if extinct

    if (successfulBirthsCount == 0 && liveCells.size > 0) {
      // console.log("Stasis detected or imminent.");
    }
}


// =============================================================================
// --- Drawing Functions (Cell drawing updated) ---
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
         drawHexagon(targetBuffer, q, r, cellCol);
    }
}

function drawLiveCellsOnBuffer(targetBuffer){ // Used for final state rendering
    if (!targetBuffer) return;
    targetBuffer.push();
    targetBuffer.translate(targetBuffer.width / 2, targetBuffer.height / 2);
    targetBuffer.clear(); // Clear previous final state
    drawLiveCells(targetBuffer);
    targetBuffer.pop();
 }

// =============================================================================
// --- UI Setup & Updates (HEAVILY MODIFIED) ---
// =============================================================================
function createUI() {
    const uiPanel = select('#ui-panel');
    uiElements = { controlGroups: {}, sliders: {}, selectors: {}, buttons: {} }; // Reset
    uiParamDisplays = {};

    // --- Helper Functions ---
    const createSliderRow = (label, paramKey, minVal, maxVal, stepVal, formatDigits = 2, groupName = null, subKey = null) => {
        let rowDiv = createDiv(label).parent(uiPanel).addClass('ui-row');
        let slider = createSlider(minVal, maxVal, subKey ? PARAMS[paramKey][subKey] : PARAMS[paramKey], stepVal)
            .parent(rowDiv)
            .input(() => {
                let val = slider.value();
                if (subKey) PARAMS[paramKey][subKey] = val;
                else PARAMS[paramKey] = val;
                select('#' + (subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay').html(nf(val, 0, formatDigits));
                if ((paramKey.includes('gaussian') || paramKey.includes('gradient') || paramKey.includes('noise') || paramKey.includes('stripe') || paramKey.includes('Band')) && !start) {
                    console.log(`Gradient/Nutrient param ${paramKey}${subKey ? '.'+subKey : ''} changed, re-init nutrients...`);
                    initGridAndNutrients(); redraw(); updateUIParamDisplays();
                }
            });
        createSpan(nf(subKey ? PARAMS[paramKey][subKey] : PARAMS[paramKey], 0, formatDigits))
            .parent(rowDiv)
            .id((subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay');
        
        uiElements.sliders[(subKey ? paramKey + '_' + subKey : paramKey)] = slider;
        uiParamDisplays[(subKey ? paramKey + '_' + subKey : paramKey)] = select('#' + (subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay');

        if (groupName) {
            if (!uiElements.controlGroups[groupName]) uiElements.controlGroups[groupName] = [];
            uiElements.controlGroups[groupName].push(rowDiv);
        }
        return rowDiv;
    };

    const createSectionHeader = (text, groupName = null) => {
        let headerP = createP(`<b>${text}</b>`).parent(uiPanel);
        if (groupName) {
            if (!uiElements.controlGroups[groupName]) uiElements.controlGroups[groupName] = [];
            uiElements.controlGroups[groupName].push(headerP);
        }
        return headerP;
    };

     const createSelectRow = (label, paramKey, options, groupName = null) => {
        let rowDiv = createDiv(label + ": ").parent(uiPanel).addClass('ui-row');
        let sel = createSelect().parent(rowDiv);
        options.forEach(opt => sel.option(typeof opt === 'string' ? opt : opt.label, typeof opt === 'string' ? opt : opt.value));
        sel.selected(PARAMS[paramKey]);
        sel.changed(() => {
            PARAMS[paramKey] = sel.value();
            console.log(`${paramKey} changed to ${sel.value()}`);
            if ((paramKey.includes('gradient') || paramKey.includes('noiseBase') || paramKey.includes('initialCondition')) && !start) {
                initGridAndNutrients(); redraw(); updateUIParamDisplays();
            } else if (paramKey === 'gradientMode') { // Special handling for mode change
                 updateUIParamDisplays(); // Show/hide relevant controls
                 if(!start) {initGridAndNutrients(); redraw();}
            }
        });
        uiElements.selectors[paramKey] = sel;
         if (groupName) {
            if (!uiElements.controlGroups[groupName]) uiElements.controlGroups[groupName] = [];
            uiElements.controlGroups[groupName].push(rowDiv);
        }
        return rowDiv;
    };


    // --- Gradient Mode ---
    createSectionHeader("Environment Type:");
    createSelectRow("Nutrient Mode", 'gradientMode',
        ['dynamic_banded_radial', 'smooth_radial', 'smooth_linear_h', 'smooth_linear_v', 'gaussian_source', 'noisy', 'stripes_h', 'stripes_v', 'uniform_g', 'uniform_h', 'uniform_mix']
    );

    // --- Dynamic Banded Radial Controls ---
    createSectionHeader("Dynamic Banded Radial Controls:", "dynamic_banded_radial");
    createSliderRow("Default Glu:", 'nutrientBandDefaultGlu', 0.0, 1.0, 0.01, 2, "dynamic_banded_radial");
    createSliderRow("Default Gal:", 'nutrientBandDefaultGal', 0.0, 1.0, 0.01, 2, "dynamic_banded_radial");
    createSliderRow("High Glu in Band:", 'nutrientBandHighGlu', 0.0, 1.0, 0.01, 2, "dynamic_banded_radial");
    createSliderRow("High Gal in Band:", 'nutrientBandHighGal', 0.0, 1.0, 0.01, 2, "dynamic_banded_radial");
    createSliderRow("Band Update Interval (cells):", 'nutrientBandChangeIntervalCells', 50, 1000, 10, 0, "dynamic_banded_radial");
    createSliderRow("Band Max Radius Factor:", 'nutrientBandMaxRadiusFactor', 1.0, 3.0, 0.1, 1, "dynamic_banded_radial");
    createSliderRow("Band Width Min (hex R):", 'nutrientBandWidthMin', 1.0, 10.0, 0.5, 1, "dynamic_banded_radial");
    createSliderRow("Band Width Max (hex R):", 'nutrientBandWidthMax', 2.0, 15.0, 0.5, 1, "dynamic_banded_radial");
    createSliderRow("Band Min Start Radius:", 'nutrientBandMinStartRadius', 0.0, 20.0, 0.5, 1, "dynamic_banded_radial");
    createSliderRow("Band Min Separation:", 'nutrientBandMinSeparation', 0.0, 10.0, 0.5, 1, "dynamic_banded_radial");


    // --- Smooth Radial/Linear Controls ---
    createSectionHeader("Smooth Radial/Linear Controls:", "smooth_radial_linear");
    createSliderRow("Center (%):", 'gradientCenter', 0.0, 1.0, 0.01, 2, "smooth_radial_linear");
    createSliderRow("Transition Width (%):", 'gradientTransitionWidth', 0.05, 1.0, 0.01, 2, "smooth_radial_linear");

    // --- Gaussian Controls ---
    createSectionHeader("Gaussian Source Controls:", "gaussian_source");
    createSliderRow("Source Q:", 'gaussianSourceQ', -maxWorldRadius, maxWorldRadius, 1, 0, "gaussian_source");
    createSliderRow("Source R:", 'gaussianSourceR', -maxWorldRadius, maxWorldRadius, 1, 0, "gaussian_source");
    createSliderRow("Width (Hex):", 'gaussianWidth', 1, maxWorldRadius/2, 0.5, 1, "gaussian_source");
    createSliderRow("Peak G:", 'gaussianPeakG', 0.0, 1.0, 0.01, 2, "gaussian_source");
    createSliderRow("Peak H:", 'gaussianPeakH', 0.0, 1.0, 0.01, 2, "gaussian_source");
    createSliderRow("Base G:", 'gaussianBaseG', 0.0, 1.0, 0.01, 2, "gaussian_source");
    createSliderRow("Base H:", 'gaussianBaseH', 0.0, 1.0, 0.01, 2, "gaussian_source");

    // --- Noisy Controls ---
    createSectionHeader("Noisy Gradient Controls:", "noisy");
    createSliderRow("Noise Scale:", 'noiseScale', 0.01, 0.5, 0.01, 2, "noisy");
    createSliderRow("Noise Inten. G:", 'noiseIntensityG', 0.0, 0.5, 0.01, 2, "noisy");
    createSliderRow("Noise Inten. H:", 'noiseIntensityH', 0.0, 0.5, 0.01, 2, "noisy");
    createSelectRow("Noise Base Grad:", 'noiseBaseGradient',
        ['smooth_radial', 'smooth_linear_h', 'smooth_linear_v', 'gaussian_source', 'uniform_g', 'uniform_h', 'uniform_mix'], "noisy"
    );

    // --- Stripe Controls ---
    createSectionHeader("Stripe Controls:", "stripes");
    createSliderRow("Stripe Width (px):", 'stripeWidth', 5, 100, 1, 0, "stripes");

    // --- Cell Genotype, Phenotype, Lag & Memory ---
    createSectionHeader("Cell Genotype, Phenotype & Lag:");
    createSelectRow("Initial Condition:", 'initialConditionType', ["mix", "all_green", "all_red"]);
    createSliderRow("Initial Green Fraction (if mix):", 'initialGenotypeMix', 0, 1, 0.01, 2);
    createSliderRow("Mut. Green -> Red:", 'mutationRateGreenToRed', 0, 0.01, 0.0001, 4);
    createSliderRow("Mut. Red -> Green:", 'mutationRateRedToGreen', 0, 0.01, 0.0001, 4);
    createSliderRow("Max Lag Time:", 'maxLagTime', 0, 2.0, 0.01, 2);
    createSliderRow("Memory Decay (gens):", 'memoryDecayGens', 0.1, 10.0, 0.1, 1);
    createSliderRow("Min Lag w/ Memory:", 'minLagTimeWithMemory', 0, 1.0, 0.01, 2);
    createSliderRow("Lag Stochasticity:", 'lagTimeStochasticity', 0, 0.5, 0.01, 2);
    createSliderRow("Time Step / Cell (dt):", 'dtPerCell', 0.001, 0.1, 0.001, 3);
    createSliderRow("Nutrient Threshold (Presence):", 'nutrientThresholdForPresence', 0.01, 0.99, 0.01, 2);
    createSliderRow("Glucose Threshold (Repression):", 'glucoseThresholdForRepression', 0.01, 0.99, 0.01, 2);


    // --- Fitness Parameters ---
    createSectionHeader("Fitness Parameters (Growth Rates):");
    const fp = 'fitnessParams'; // Alias for brevity
    createSliderRow("Green, Gal_ON, on Gal:", fp, 0, 3.0, 0.05, 2, null, 'fitness_green_on_gal');
    createSliderRow("Red, Gal_ON, on Gal:", fp, 0, 3.0, 0.05, 2, null, 'fitness_red_on_gal');
    createSliderRow("Gal_ON, no Gal (cost factor):", fp, 0, 1.0, 0.01, 2, null, 'fitness_on_no_gal_cost_factor');
    createSliderRow("Green, Gal_OFF, on Glu:", fp, 0, 3.0, 0.05, 2, null, 'fitness_green_off_glu');
    createSliderRow("Red, Gal_OFF, on Glu:", fp, 0, 3.0, 0.05, 2, null, 'fitness_red_off_glu');
    createSliderRow("Gal_OFF, on Gal (inducing):", fp, 0, 0.5, 0.01, 2, null, 'fitness_off_on_gal');
    createSliderRow("Fitness during Lag:", fp, 0, 0.5, 0.01, 2, null, 'fitness_during_lag');


    // --- Colonization Parameters ---
    createSectionHeader("Colonization Parameters:");
    createSliderRow("Min Viability Fitness:", 'minColonizationFitness', 0, 0.5, 0.01, 2);
    createSliderRow("Pressure Factor:", 'colonizationPressureFactor', 0.1, 2.0, 0.1, 1);
    createSliderRow("Max Colonize Prob:", 'maxColonizationProbability', 0.5, 1.0, 0.01, 2);
    createSliderRow("Initial Patch Radius:", 'initialPatchRadius', 0, 10, 1, 0);


    // --- Controls Section ---
    createSectionHeader("Controls:");
    let controlsDiv = createDiv().parent(uiPanel).id('controls-container');
    uiElements.buttons.startBtn = createButton("Start / Pause (Space)").parent(controlsDiv).mousePressed(toggleSimulation);
    uiElements.buttons.resetBtn = createButton("Reset Sim").parent(controlsDiv).mousePressed(() => {
        console.log("Resetting simulation..."); initGridAndNutrients(); redraw();
    });
    uiElements.buttons.gradientOverlayBtn = createButton("Show Gradient").parent(controlsDiv).mousePressed(toggleGradientOverlay);
    uiElements.gradientOverlayBtn = uiElements.buttons.gradientOverlayBtn; // for compatibility with old reference


    // --- Status Display Area ---
    let statusDiv = createDiv().parent(uiPanel).id('status-section');
    createP("<b>Status:</b>").parent(statusDiv);
    createSpan('Step: 0').parent(statusDiv).id('stepValueDisplay');
    createSpan('Cell Count: 0').parent(statusDiv).id('cellCountValueDisplay');
    createSpan('Max Radius: 0').parent(statusDiv).id('maxRadiusValueDisplay');
    createSpan('Status: Initialized').parent(statusDiv).id('statusValueDisplay');

    updateUIParamDisplays(); // Initial setup of visibility & values
}


function toggleSimulation() {
    start = !start;
    if (start) { loop(); console.log("Sim Started."); }
    else { noLoop(); console.log("Sim Paused."); redraw(); } // Redraw once when paused
    updateDynamicUIText();
}

function toggleGradientOverlay() {
    showGradientOverlay = !showGradientOverlay;
    uiElements.gradientOverlayBtn.html(showGradientOverlay ? 'Show Cells' : 'Show Gradient');
    if (!start) redraw(); // Redraw to show change if paused
}

function updateUIParamDisplays() {
    // Update slider values displayed next to them
    for (const key in uiParamDisplays) {
        let valueToDisplay;
        let pKey = key;
        let subKey = null;
        if (key.includes('_')) { // handles fitnessParams.sub_key
            const parts = key.split('_');
            pKey = parts[0];
            subKey = parts.slice(1).join('_'); // Rejoin if subkey has underscores
             if (PARAMS[pKey] && typeof PARAMS[pKey] === 'object' && PARAMS[pKey].hasOwnProperty(subKey)) {
                valueToDisplay = PARAMS[pKey][subKey];
            } else {continue;} // Should not happen if UI setup correctly
        } else if (PARAMS.hasOwnProperty(key)) {
            valueToDisplay = PARAMS[key];
        } else {
            continue; // Param not found
        }

        if (uiParamDisplays[key]) {
            // Infer formatDigits from slider step or default
            let digits = 2;
            if (uiElements.sliders[key]) {
                let step = uiElements.sliders[key].elt.step;
                if (step.includes('.')) digits = step.split('.')[1].length;
                else if (parseFloat(step) === 1) digits = 0;
            }
             // Override for specific known params
            if (key === 'mutationRateGreenToRed' || key === 'mutationRateRedToGreen') digits = 4;
            else if (key === 'dtPerCell') digits = 3;
            else if (key === 'nutrientBandChangeIntervalCells' || key === 'initialPatchRadius') digits = 0;
            else if (key === 'nutrientBandMaxRadiusFactor' || key === 'memoryDecayGens' ||
                     key === 'colonizationPressureFactor' || key === 'gaussianWidth' ||
                     key === 'nutrientBandWidthMin' || key === 'nutrientBandWidthMax' ||
                     key === 'nutrientBandMinStartRadius' || key === 'nutrientBandMinSeparation') digits = 1;


            if (uiParamDisplays[key].html() !== nf(valueToDisplay, 0, digits)) {
                 uiParamDisplays[key].html(nf(valueToDisplay, 0, digits));
            }
        }
    }

    // Update visibility of gradient controls
    const currentMode = PARAMS.gradientMode;
    const controlGroupMapping = {
        "dynamic_banded_radial": "dynamic_banded_radial",
        "smooth_radial": "smooth_radial_linear",
        "smooth_linear_h": "smooth_radial_linear",
        "smooth_linear_v": "smooth_radial_linear",
        "gaussian_source": "gaussian_source",
        "noisy": "noisy",
        "stripes_h": "stripes",
        "stripes_v": "stripes"
        // uniform modes have no specific controls
    };

    for (const groupName in uiElements.controlGroups) {
        let showGroup = (controlGroupMapping[currentMode] === groupName);
        uiElements.controlGroups[groupName].forEach(element => {
            if (element) {
                if (showGroup) element.show(); else element.hide();
            }
        });
    }
}


function updateDynamicUIText() {
    select('#stepValueDisplay')?.html(`Step: ${simulationSteps}`);
    select('#cellCountValueDisplay')?.html(`Cell Count: ${liveCells.size}`);
    select('#maxRadiusValueDisplay')?.html(`Max Radius: ${currentColonyMaxRadius.toFixed(1)}`);
    let statusSpan = select('#statusValueDisplay');

    if (statusSpan) {
        let statusMsg = "Initialized"; let statusColor = "#6c757d"; // Default gray
        if (start) { statusMsg = "RUNNING"; statusColor = "#28a745"; } // Green
        else if (simulationSteps > 0) {
            if (liveCells.size === 0) { statusMsg = "EXTINCTION"; statusColor = "#dc3545"; } // Red
            else if (lastStepBirths === 0 && finalStateRendered) { statusMsg = "STASIS - Stopped"; statusColor = "#dc3545"; } // Red
            else if (lastStepBirths === 0) { statusMsg = "PAUSED (Stasis likely)"; statusColor = "#ffc107"; } // Yellow
            else { statusMsg = "PAUSED"; statusColor = "#ffc107"; } // Yellow
        } else { statusMsg = "Ready / PAUSED"; statusColor = "#6c757d"; } // Gray
        statusSpan.html(`Status: ${statusMsg}`).style('color', statusColor);
    }
}

// =============================================================================
// --- User Interaction (Planting updated) ---
// =============================================================================
function mousePressed() {
    if (mouseX >= 0 && mouseX <= width && mouseY >= 0 && mouseY <= height) {
        if (mouseX > canvasWidth) return; // Ignore clicks in UI panel

        let clickX_in_grid_space = mouseX - gridOriginX;
        let clickY_in_grid_space = mouseY - gridOriginY;
        const hexCoords = pixelToHex(clickX_in_grid_space, clickY_in_grid_space);

        if (axialDistance(hexCoords.q, hexCoords.r, 0, 0) <= maxWorldRadius + 2) {
            plantPatch(hexCoords.q, hexCoords.r, PARAMS.initialPatchRadius);
            if (!start) { redraw(); }
        } else {
             console.log("Click outside defined world radius.");
        }
    }
}

function plantPatch(centerQ, centerR, radius) {
    let count = 0; let newlyPlantedCoords = [];
    const { initialConditionType, initialGenotypeMix, memoryDecayGens } = PARAMS;

    for (let q = -radius; q <= radius; q++) {
        let r1 = Math.max(-radius, -q - radius);
        let r2 = Math.min(radius, -q + radius);
        for (let r = r1; r <= r2; r++) {
            let targetQ = centerQ + q;
            let targetR = centerR + r;
            const coordStr = `${targetQ},${targetR}`;
            if (glucose.has(coordStr) && !liveCells.has(coordStr)) { // Check if hex exists in nutrient grid and is empty
                let genotype_val = 0; // Default Red
                if (initialConditionType === "all_green") genotype_val = 1;
                else if (initialConditionType === "mix" && random() < initialGenotypeMix) genotype_val = 1;

                const newCell = {
                    genotype: genotype_val,
                    galPhenotype: "Gal_OFF", // Start OFF
                    gensSinceGalOn: memoryDecayGens * 5, // effectively no memory
                    lagTimeRemaining: 0.0
                };
                liveCells.set(coordStr, newCell);
                newlyPlantedCoords.push(coordStr);
                count++;
                currentColonyMaxRadius = max(currentColonyMaxRadius, axialDistance(targetQ, targetR, 0,0));
            }
        }
    }

    if (count > 0) {
        // Update frontier based on newly planted cells
        let cellsToConsiderForFrontier = new Set(newlyPlantedCoords);
        // Also consider neighbors of newly planted cells if they were already live
        for (const coordStr of newlyPlantedCoords) {
            const [q_p, r_p] = coordStr.split(',').map(Number);
            for (const dir of axialDirections) {
                const nq = q_p + dir.q; const nr = r_p + dir.r;
                const neighborCoordStr = `${nq},${nr}`;
                if (liveCells.has(neighborCoordStr) && !newlyPlantedCoords.includes(neighborCoordStr)) {
                    cellsToConsiderForFrontier.add(neighborCoordStr);
                }
            }
        }

        for (const coordStrToEvaluate of cellsToConsiderForFrontier) {
            const [q_e, r_e] = coordStrToEvaluate.split(',').map(Number);
            let hasEmptyNeighbor = false;
            for (const dir of axialDirections) {
                const nq = q_e + dir.q; const nr = r_e + dir.r;
                const neighborCoordStr = `${nq},${nr}`;
                if (glucose.has(neighborCoordStr) && !liveCells.has(neighborCoordStr)) {
                    hasEmptyNeighbor = true; break;
                }
            }
            if (hasEmptyNeighbor) frontierSet.add(coordStrToEvaluate);
            else frontierSet.delete(coordStrToEvaluate); // Should not happen for newly planted but good for existing neighbors
        }

        console.log(`Planted ${count} cells near (q:${centerQ}, r:${centerR}). Initial type: ${initialConditionType}. Frontier updated.`);
        lastStepBirths = count; // To avoid immediate stasis if paused
        if (!start && finalStateRendered) {
            console.log("Planting invalidated pre-rendered final state.");
            finalStateRendered = false;
        }
        if (!start) { updateDynamicUIText(); redraw(); }
    } else {
        console.log(`Could not plant near (q:${centerQ}, r:${centerR}) - occupied or outside range.`);
    }
}


function keyPressed() {
  if (key === ' ') { toggleSimulation(); return false; } // Prevent browser scroll
}

// =============================================================================
