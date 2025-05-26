// --- Core Microbial Strategy Simulation (v32.2 - Perf/Mem Opt & Logic Review) ---
// Model based on "2D Radial Microbial Range Expansion in Fluctuating Environments" description
// Focus: Performance, Memory, Adherence to original model logic.
// ---

// --- Simulation Settings ---
let worldSize = 120; // Half-width of the square world (e.g., -worldSize to +worldSize)
let cellSize = 4;    // Pixel size of each cell (smaller for more cells)
let dt = 0.1;        // Simulation time step
const epsilon = 1e-9; // Small number for float comparisons
let drawInterval = 1; // Draw every frame

// --- MODEL PARAMETERS (Defaults, adjustable via UI) ---
let PARAMS = {
  // --- Environment ---
  W_band: 15, // Width of each nutrient band (in lattice units)

  // --- Cell General ---
  initial_cluster_radius: 3,
  default_strategy_type: "Responsive",
  lambda_L_factor: 0.7, // Growth rate on Galactose = lambda_L_factor * inherent_growth_rate_G
  max_possible_growth_G: 1.0, // For normalization in trade-offs/lag calcs
  initial_inherent_growth_rate_G: 0.8,

  // --- Trade-off: Growth Rate G vs. Lag Time G->L ---
  trade_off_base_lag: 0.5,
  trade_off_max_additional_lag: 3.0,
  trade_off_exponent_k: 2.0,

  // --- DK-Inspired Lag Modification for New Cells ---
  lag_mod_factor: 0.5, // Positive: better G-parents increase daughter's G->L lag

  // --- Strategy: Constitutive ---
  constitutive_growth_G: 0.6,
  constitutive_growth_L: 0.5,

  // --- Strategy: BetHedging ---
  bet_hedging_p_switch_to_L_if_G_parent: 0.1,
  bet_hedging_p_switch_to_G_if_L_parent: 0.1,

  // --- Strategy: Memory ---
  memory_M_gen_memory_window: 5.0,
  memory_T_lag_GL_multiplier: 0.2,
};

// --- Simulation State ---
let nutrientGrid = new Map(); // Stores "Glucose" or "Galactose" for each (x,y) string key "x,y"
let liveCells = new Map();    // Stores cell objects, key is "x,y" string
let nextCellID = 0;
let simulationTime = 0;
let simulationSteps = 0;
let maxColonyRadius = 0;
let previousMaxColonyRadiusForVelocity = 0;
let timeAtPreviousMaxRadiusForVelocity = 0;
let avgRadialVelocity = 0;
let lastStepBirths = 1;       // To detect stasis
let start = false;
let showNutrientOverlay = false;

// --- Graphics & Grid ---
let nutrientBuffer, cellBuffer, finalStateBuffer;
let finalStateRendered = false;
let canvasWidth, canvasHeight, bufferWidth, bufferHeight, gridOriginX, gridOriginY;

const moore_deltas = [
  { dx: -1, dy: -1 }, { dx: 0, dy: -1 }, { dx: 1, dy: -1 },
  { dx: -1, dy: 0 },  /* {dx:0, dy:0} */ { dx: 1, dy: 0 },
  { dx: -1, dy: 1 },  { dx: 0, dy: 1 },  { dx: 1, dy: 1 },
];

// --- Color Definitions ---
let colorGlucoseBand, colorGalactoseBand;
let colorGSpecialist, colorLSpecialist, colorSwitchingGL;

// --- UI Elements ---
let uiElements = { controlGroups: {}, sliders: {}, selectors: {}, buttons: {} };
let uiParamDisplays = {};
let paramMapping = {}; // For robust UI to PARAMS mapping


// =============================================================================
// --- Math & Grid Helpers ---
// =============================================================================
function euclideanDistance(x1, y1, x2, y2) {
  return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
}

function pixelToGridCoords(pixelX, pixelY) {
  let worldX = pixelX - gridOriginX; // gridOriginX is buffer center, so pixelX is canvas coord
  let worldY = pixelY - gridOriginY; // This gives coords relative to grid (0,0) pixel space
  return { 
    x: Math.floor((worldX + cellSize / 2) / cellSize), 
    y: Math.floor((worldY + cellSize / 2) / cellSize) 
  };
}

function gridToPixelCoords(gridX, gridY) {
  // Top-left pixel of the cell, relative to the (0,0) of the grid in pixel space
  return { 
    x: gridX * cellSize - cellSize / 2, 
    y: gridY * cellSize - cellSize / 2 
  };
}

// =============================================================================
// --- Drawing Helpers ---
// =============================================================================
function calculateGridDimensions() {
  bufferWidth = (worldSize * 2 + 1) * cellSize;
  bufferHeight = (worldSize * 2 + 1) * cellSize;
  canvasWidth = bufferWidth + 40; // some padding
  canvasHeight = bufferHeight + 40; // Increased for status text below canvas
  gridOriginX = bufferWidth / 2; // Relative to buffer's own top-left (i.e., buffer center)
  gridOriginY = bufferHeight / 2; // Relative to buffer's own top-left
}

function drawSquareCell(buffer, gridX, gridY, fillColor, strokeColor = null, weight = 0.5) {
  const p = gridToPixelCoords(gridX, gridY); // Gets top-left relative to grid (0,0) in pixel space
  buffer.fill(fillColor);
  if (strokeColor) {
    buffer.stroke(strokeColor);
    buffer.strokeWeight(weight);
  } else {
    buffer.noStroke();
  }
  buffer.rect(p.x , p.y , cellSize, cellSize);
}

// =============================================================================
// --- P5.js Core Functions ---
// =============================================================================
function setup() {
  calculateGridDimensions();
  nutrientBuffer = createGraphics(bufferWidth, bufferHeight);
  cellBuffer = createGraphics(bufferWidth, bufferHeight);
  finalStateBuffer = createGraphics(bufferWidth, bufferHeight);
  finalStateBuffer.pixelDensity(1);

  colorGlucoseBand = color(255, 230, 230, 180);
  colorGalactoseBand = color(230, 230, 255, 180);
  colorGSpecialist = color(0, 200, 0);
  colorLSpecialist = color(0, 0, 200);
  colorSwitchingGL = color(255, 165, 0);

  let totalCanvasHeightUI = Math.max(canvasHeight, 750); // For UI panel height
  let mainCanvas = createCanvas(canvasWidth, totalCanvasHeightUI);
  let canvasContainer = select('#canvas-container');
  if (canvasContainer) mainCanvas.parent('canvas-container');

  createUI();
  initSimulation();

  background(240);
  push();
  // The main canvas origin for drawing buffers is top-left for image()
  // The buffers themselves handle their internal translations via gridOriginX/Y
  translate((canvasWidth - bufferWidth)/2, (canvasHeight - bufferHeight)/2);
  image(nutrientBuffer, 0, 0);
  image(cellBuffer, 0, 0);
  pop();

  updateUIParamDisplays();
  updateDynamicUIText();
  noLoop();
}

function draw() {
  if (start) {
    stepSimulation();
  }

  if (frameCount % drawInterval === 0 || !start || showNutrientOverlay) {
    background(240); // Clear main canvas
    push();
    translate((canvasWidth - bufferWidth)/2, (canvasHeight - bufferHeight)/2); // Center buffers on canvas

    if (showNutrientOverlay) {
      image(nutrientBuffer, 0, 0);
    } else {
      image(nutrientBuffer, 0, 0); // Always show nutrients underneath
      if (!start && finalStateRendered) {
        image(finalStateBuffer, 0, 0);
      } else {
        // cellBuffer is updated in stepSimulation if running, or initSimulation
        image(cellBuffer, 0, 0);
      }
    }
    pop();
  }
  updateDynamicUIText();
}

// =============================================================================
// --- Nutrient Environment Logic ---
// =============================================================================
function getNutrientAt(x, y) { // Grid coordinates
  const dist = euclideanDistance(x, y, 0, 0);
  const band_index = Math.floor(dist / PARAMS.W_band);
  return (band_index % 2 === 0) ? "Glucose" : "Galactose";
}

function populateNutrientGridAndDrawBuffer() {
  nutrientGrid.clear();
  nutrientBuffer.push();
  nutrientBuffer.translate(gridOriginX, gridOriginY); // Translate to grid (0,0) for drawing cells
  nutrientBuffer.background(220); // Default background for out-of-bounds

  for (let gx = -worldSize; gx <= worldSize; gx++) {
    for (let gy = -worldSize; gy <= worldSize; gy++) {
      const nutrientType = getNutrientAt(gx, gy);
      nutrientGrid.set(`${gx},${gy}`, nutrientType);
      const fillColor = (nutrientType === "Glucose") ? colorGlucoseBand : colorGalactoseBand;
      drawSquareCell(nutrientBuffer, gx, gy, fillColor, null); // No stroke for nutrient background
    }
  }
  nutrientBuffer.pop();
}


// =============================================================================
// --- Cell Agent Logic ---
// =============================================================================
function trade_off_lag_vs_growth(inherent_growth_rate_G) {
  const { trade_off_base_lag, trade_off_max_additional_lag, trade_off_exponent_k, max_possible_growth_G } = PARAMS;
  if (inherent_growth_rate_G <= epsilon) return trade_off_base_lag + trade_off_max_additional_lag;
  const normalized_gr_G = Math.min(inherent_growth_rate_G, max_possible_growth_G) / max_possible_growth_G;
  let lag = trade_off_base_lag + Math.pow(normalized_gr_G, trade_off_exponent_k) * trade_off_max_additional_lag;
  return Math.max(0, lag); // Ensure lag is not negative
}

function calculate_modified_lag(inherent_T_lag_GL_daughter, avg_inherent_G_growth_parents) {
  const { lag_mod_factor, max_possible_growth_G } = PARAMS;
  if (avg_inherent_G_growth_parents < epsilon) return inherent_T_lag_GL_daughter;
  const normalized_avg_fitness = Math.min(avg_inherent_G_growth_parents, max_possible_growth_G) / max_possible_growth_G;
  let modifiedLag = inherent_T_lag_GL_daughter * (1 + lag_mod_factor * normalized_avg_fitness);
  return Math.max(0, modifiedLag); // Ensure lag is not negative
}

function createCell(x, y, strategy_type, inherent_growth_rate_G_val, parentCell = null) {
  nextCellID++;
  const inherent_T_lag_GL_val = trade_off_lag_vs_growth(inherent_growth_rate_G_val);
  let newCell = {
    id: nextCellID, x: x, y: y,
    phenotype: "G_specialist", // Default, will be adjusted immediately after creation
    strategy_type: strategy_type,
    inherent_growth_rate_G: inherent_growth_rate_G_val,
    inherent_T_lag_GL: inherent_T_lag_GL_val,
    current_growth_rate: 0.0,
    remaining_lag_time: 0.0,
    // Inherit parent's memory state and advance it by dt, or assume no recent memory if no parent
    memory_state: parentCell ? parentCell.memory_state + dt : PARAMS.memory_M_gen_memory_window + 1,
  };
  return newCell;
}


// =============================================================================
// --- Initialization & Reset ---
// =============================================================================
function initSimulation() {
  start = false; loop(); // enable draw loop briefly for reset functions to take effect
  lastStepBirths = 1; finalStateRendered = false; showNutrientOverlay = false;
  if (uiElements.buttons?.nutrientOverlayBtn) uiElements.buttons.nutrientOverlayBtn.html('Show Nutrients');
  if (finalStateBuffer) finalStateBuffer.clear();
  cellBuffer.clear(); // Clear cell drawing buffer

  simulationTime = 0; simulationSteps = 0; maxColonyRadius = 0; nextCellID = 0;
  previousMaxColonyRadiusForVelocity = 0; timeAtPreviousMaxRadiusForVelocity = 0; avgRadialVelocity = 0;

  liveCells.clear();
  populateNutrientGridAndDrawBuffer(); // Sets up nutrientGrid and nutrientBuffer

  const { initial_cluster_radius, default_strategy_type, initial_inherent_growth_rate_G } = PARAMS;
  let plantedCount = 0;
  for (let gx = -initial_cluster_radius; gx <= initial_cluster_radius; gx++) {
    for (let gy = -initial_cluster_radius; gy <= initial_cluster_radius; gy++) {
      if (euclideanDistance(gx, gy, 0, 0) <= initial_cluster_radius) {
        if (nutrientGrid.has(`${gx},${gy}`)) { // Ensure it's within world bounds
          const coordStr = `${gx},${gy}`;
          if (!liveCells.has(coordStr)) {
            let cell = createCell(gx, gy, default_strategy_type, initial_inherent_growth_rate_G);
            const local_nutrient = nutrientGrid.get(coordStr);
            
            if (local_nutrient === "Glucose") {
                cell.phenotype = "G_specialist";
                cell.current_growth_rate = cell.inherent_growth_rate_G;
                if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
                    cell.memory_state = PARAMS.memory_M_gen_memory_window + 1; // No recent Gal memory
                }
            } else { // Started on Galactose
                cell.phenotype = "Switching_GL";
                cell.remaining_lag_time = cell.inherent_T_lag_GL;
                cell.current_growth_rate = 0.0;
                if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
                    cell.memory_state = 0; // Starting on Gal, so fresh memory of it
                }
            }
            liveCells.set(coordStr, cell);
            maxColonyRadius = Math.max(maxColonyRadius, euclideanDistance(gx, gy, 0, 0));
            plantedCount++;
          }
        }
      }
    }
  }
  previousMaxColonyRadiusForVelocity = maxColonyRadius; // Initialize for velocity calc

  console.log(`Initialized ${plantedCount} cells. Strategy: ${default_strategy_type}. Max Radius: ${maxColonyRadius.toFixed(2)}`);
  drawLiveCellsOnBuffer(cellBuffer); // Initial draw
  noLoop(); // Start paused
  updateDynamicUIText();
}

// =============================================================================
// --- Simulation Step Logic ---
// =============================================================================
function stepSimulation() {
  if (liveCells.size === 0 || (lastStepBirths === 0 && simulationSteps > 10)) { // Added step condition for stasis
    if (start) {
      let reason = liveCells.size === 0 ? "Extinction" : "Stasis";
      console.log(`${reason} @ step ${simulationSteps}, t=${simulationTime.toFixed(2)}. Stopping.`);
      if (liveCells.size > 0 && !finalStateRendered) {
        drawLiveCellsOnBuffer(finalStateBuffer); finalStateRendered = true;
      }
      start = false; noLoop(); updateDynamicUIText();
    }
    return;
  }

  simulationSteps++;
  simulationTime += dt;
  let currentBirthsThisStep = 0;

  // --- 1. Cell State Update (for existing cells) ---
  for (const cell of liveCells.values()) { // Iterate directly over map values
    const local_nutrient = nutrientGrid.get(`${cell.x},${cell.y}`);
    let switchedToLagThisStep = false;

    // A. Lag countdown
    if (cell.phenotype === "Switching_GL") {
      cell.remaining_lag_time -= dt;
      if (cell.remaining_lag_time <= epsilon) {
        cell.remaining_lag_time = 0;
        cell.phenotype = "L_specialist";
        if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
             cell.memory_state = 0; // Switched to L, reset memory timer
        }
      }
    }
    // B. Strategy-specific phenotype switching (if not lagging)
    else {
      if (cell.strategy_type === "Responsive" || cell.strategy_type === "BetHedging" || cell.strategy_type === "Memory") {
        if (cell.phenotype === "G_specialist" && local_nutrient === "Galactose") {
          cell.phenotype = "Switching_GL";
          switchedToLagThisStep = true;
          if ((cell.strategy_type === "Memory" || cell.strategy_type === "Responsive") && PARAMS.memory_M_gen_memory_window > 0) {
            if (cell.memory_state < PARAMS.memory_M_gen_memory_window) {
              cell.remaining_lag_time = cell.inherent_T_lag_GL * PARAMS.memory_T_lag_GL_multiplier;
            } else {
              cell.remaining_lag_time = cell.inherent_T_lag_GL;
            }
          } else { // BetHedging or Responsive without active memory effect
            cell.remaining_lag_time = cell.inherent_T_lag_GL;
          }
        } else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") {
          cell.phenotype = "G_specialist"; // No L->G lag
        }
      }
      if (switchedToLagThisStep && (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0))) {
            cell.memory_state = 0; // Just started switch to L due to Gal, memory is "fresh"
       }
    }

    // C. Update current_growth_rate
    cell.current_growth_rate = 0.0;
    if (cell.phenotype === "Switching_GL") { /* Growth is 0 */ }
    else if (cell.strategy_type === "Constitutive") {
      cell.current_growth_rate = (local_nutrient === "Glucose") ? PARAMS.constitutive_growth_G : PARAMS.constitutive_growth_L;
    } else { // For Responsive, Memory, BetHedging (once phenotype set)
      if (cell.phenotype === "G_specialist" && local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G;
      else if (cell.phenotype === "L_specialist" && local_nutrient === "Galactose") cell.current_growth_rate = PARAMS.lambda_L_factor * cell.inherent_growth_rate_G;
      else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G; // L on G, grows fine (pheno would have switched to G already)
      else if (cell.phenotype === "G_specialist" && local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G; // Explicitly for G_specialist on G
    }

    // D. Update memory_state
    if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
      if (!((cell.phenotype === "L_specialist" && local_nutrient === "Galactose") || 
            (cell.phenotype === "Switching_GL" && local_nutrient === "Galactose"))) {
        cell.memory_state += dt;
      }
      // Reset of memory_state to 0 is handled when phenotype becomes L_specialist on Gal or starts Switching_GL on Gal.
    }
  }

  // --- 2. New Cell Generation (at Frontier) ---
  // Pass 1: Identify all unique empty frontier slots.
  let emptyFrontierSlotCoords = new Set();
  for (const cell of liveCells.values()) {
      for (const delta of moore_deltas) {
          const nx = cell.x + delta.dx;
          const ny = cell.y + delta.dy;
          // Bounds check for the potential empty slot
          if (nx < -worldSize || nx > worldSize || ny < -worldSize || ny > worldSize) continue;
          
          const nCoordStr = `${nx},${ny}`;
          if (!liveCells.has(nCoordStr)) { // Is an empty slot
              emptyFrontierSlotCoords.add(nCoordStr);
          }
      }
  }

  // Pass 2: For each unique empty slot, find its parental neighbors and their properties.
  let frontierSlotsData = new Map(); // Key: "x,y" of empty slot, Value: { parentalNeighbors: [cellObj,...], sumParentGrowthRate: Z, sumParentInherentGGrowth: X, parentCount: Y }
  for (const slotCoordStr of emptyFrontierSlotCoords) {
      let slotEntry = { parentalNeighbors: [], sumParentGrowthRate: 0, sumParentInherentGGrowth: 0, parentCount:0 };
      const [sx, sy] = slotCoordStr.split(',').map(Number);

      for (const delta of moore_deltas) { // Check neighbors of the *slot*
          const px = sx + delta.dx;
          const py = sy + delta.dy;
          // Bounds check for parent position (px, py)
          if (px < -worldSize || px > worldSize || py < -worldSize || py > worldSize) continue;

          const pCoordStr = `${px},${py}`;
          if (liveCells.has(pCoordStr)) {
              const parentCell = liveCells.get(pCoordStr);
              slotEntry.parentalNeighbors.push(parentCell);
              slotEntry.sumParentGrowthRate += parentCell.current_growth_rate;
              slotEntry.sumParentInherentGGrowth += parentCell.inherent_growth_rate_G;
              slotEntry.parentCount++;
          }
      }
      if (slotEntry.parentCount > 0) { // Only consider slots that actually have parents
          frontierSlotsData.set(slotCoordStr, slotEntry);
      }
  }

  let newCellsToAdd = []; // Store {coordStr: "x,y", cell: newCellObject}
  for (const [slotCoordStr, slotData] of frontierSlotsData) { // Iterate over slots with parental data
    let chosenParent = null;
    if (slotData.sumParentGrowthRate > epsilon) {
      let randVal = random() * slotData.sumParentGrowthRate;
      let cumulativeGrowth = 0;
      for (const p of slotData.parentalNeighbors) {
        cumulativeGrowth += p.current_growth_rate;
        if (randVal <= cumulativeGrowth) { chosenParent = p; break; }
      }
      if (!chosenParent) chosenParent = slotData.parentalNeighbors[slotData.parentalNeighbors.length -1]; // Fallback if float error
    } else { // All potential parents have 0 growth rate
      chosenParent = random(slotData.parentalNeighbors);
    }
    if (!chosenParent) continue;

    const [nx, ny] = slotCoordStr.split(',').map(Number);
    let daughterCell = createCell(nx, ny, chosenParent.strategy_type, chosenParent.inherent_growth_rate_G, chosenParent);
    const local_nutrient_at_slot = nutrientGrid.get(slotCoordStr);
    daughterCell.phenotype = chosenParent.phenotype; // Initial inheritance

    if (daughterCell.strategy_type === "BetHedging") {
        if (chosenParent.phenotype === "G_specialist" && random() < PARAMS.bet_hedging_p_switch_to_L_if_G_parent) daughterCell.phenotype = "L_specialist";
        else if (chosenParent.phenotype === "L_specialist" && random() < PARAMS.bet_hedging_p_switch_to_G_if_L_parent) daughterCell.phenotype = "G_specialist";
    }
    
    let daughterSwitchedToLag = false;
    if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Galactose") {
      daughterCell.phenotype = "Switching_GL";
      daughterSwitchedToLag = true;
      let avg_inherent_G_parents = (slotData.parentCount > 0) ? slotData.sumParentInherentGGrowth / slotData.parentCount : 0;
      daughterCell.remaining_lag_time = calculate_modified_lag(daughterCell.inherent_T_lag_GL, avg_inherent_G_parents);

      if ((daughterCell.strategy_type === "Memory" || (daughterCell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) &&
          chosenParent.memory_state < PARAMS.memory_M_gen_memory_window) {
        let memoryLag = daughterCell.inherent_T_lag_GL * PARAMS.memory_T_lag_GL_multiplier;
        daughterCell.remaining_lag_time = Math.min(daughterCell.remaining_lag_time, memoryLag);
      }
    } else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Glucose") {
      daughterCell.phenotype = "G_specialist"; daughterCell.remaining_lag_time = 0;
    } else { // No switch needed or L->G (which is instant) or already correct phenotype
      daughterCell.remaining_lag_time = 0;
    }

    daughterCell.current_growth_rate = 0.0; // Default
    if (daughterCell.phenotype === "Switching_GL") { /* Growth 0 */ }
    else if (daughterCell.strategy_type === "Constitutive") daughterCell.current_growth_rate = (local_nutrient_at_slot === "Glucose") ? PARAMS.constitutive_growth_G : PARAMS.constitutive_growth_L;
    else {
      if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Glucose") daughterCell.current_growth_rate = daughterCell.inherent_growth_rate_G;
      else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Galactose") daughterCell.current_growth_rate = PARAMS.lambda_L_factor * daughterCell.inherent_growth_rate_G;
      // G_specialist on Gal (will switch / is switching, growth 0)
      // L_specialist on Glu (will switch to G_specialist, grows as G) - this case covered by phenotype switch then G growth
    }
    
    if (daughterCell.strategy_type === "Memory" || (daughterCell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
        if ((daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Galactose") || daughterSwitchedToLag) {
            daughterCell.memory_state = 0; // Born onto Gal (or switching to it), fresh memory
        }
        // else, memory_state inherited from parent and advanced in createCell
    }

    newCellsToAdd.push({ coordStr: slotCoordStr, cell: daughterCell });
    currentBirthsThisStep++;
    maxColonyRadius = Math.max(maxColonyRadius, euclideanDistance(nx, ny, 0, 0));
  }

  newCellsToAdd.forEach(item => liveCells.set(item.coordStr, item.cell));
  lastStepBirths = currentBirthsThisStep;

  if (simulationTime > timeAtPreviousMaxRadiusForVelocity + dt * 10 && maxColonyRadius > previousMaxColonyRadiusForVelocity + 0.1) {
      const dr = maxColonyRadius - previousMaxColonyRadiusForVelocity;
      const dt_vel = simulationTime - timeAtPreviousMaxRadiusForVelocity;
      if (dt_vel > epsilon) avgRadialVelocity = dr / dt_vel;
      previousMaxColonyRadiusForVelocity = maxColonyRadius;
      timeAtPreviousMaxRadiusForVelocity = simulationTime;
  }
  
  if (start) drawLiveCellsOnBuffer(cellBuffer); // Update cell buffer if running
}

// =============================================================================
// --- Drawing Functions ---
// =============================================================================
function drawLiveCellsOnBuffer(targetBuffer) {
  if (!targetBuffer) return;
  targetBuffer.clear(); // Clear previous frame on this buffer
  targetBuffer.push();
  targetBuffer.translate(gridOriginX, gridOriginY); // Center drawing

  for (const [coordStr, cell] of liveCells.entries()) {
    let cellCol;
    if (cell.phenotype === "Switching_GL") cellCol = colorSwitchingGL;
    else if (cell.phenotype === "G_specialist") cellCol = colorGSpecialist;
    else cellCol = colorLSpecialist; // L_specialist
    drawSquareCell(targetBuffer, cell.x, cell.y, cellCol, color(50), 0.2);
  }
  targetBuffer.pop();
}

// =============================================================================
// --- UI Setup & Updates ---
// =============================================================================
function createUI() {
    const uiPanel = select('#ui-panel'); uiPanel.html('');
    uiElements = { controlGroups: {}, sliders: {}, selectors: {}, buttons: {} };
    uiParamDisplays = {}; paramMapping = {}; // Reset mappings

    // Helper to map a UI key (string) to a path in PARAMS (array of strings)
    const addParamToMapping = (uiKey, ...paramPath) => {
        paramMapping[uiKey] = paramPath;
    };

    // Helper to get value from PARAMS using mapped path
    const getParamValueFromPath = (uiKey) => {
        const path = paramMapping[uiKey];
        if (!path || path.length === 0) return undefined;
        let value = PARAMS;
        for (const key of path) {
            if (value && typeof value === 'object' && value.hasOwnProperty(key)) {
                value = value[key];
            } else {
                // Fallback for flat PARAMS structure if path is just one element
                if (path.length === 1 && PARAMS.hasOwnProperty(key)) return PARAMS[key];
                return undefined; // Path is invalid or property doesn't exist
            }
        }
        return value;
    };
    
    // Helper to set value in PARAMS using mapped path
    const setParamValueByPath = (uiKey, value) => {
        const path = paramMapping[uiKey];
        if (!path || path.length === 0) return;
        let obj = PARAMS;
        for (let i = 0; i < path.length - 1; i++) {
            if (obj && typeof obj === 'object' && obj.hasOwnProperty(path[i])) {
                obj = obj[path[i]];
            } else { return; } // Path is invalid
        }
        if (obj && typeof obj === 'object' && obj.hasOwnProperty(path[path.length - 1])) {
            obj[path[path.length - 1]] = value;
        } else if (path.length === 1 && PARAMS.hasOwnProperty(path[0])) { // Fallback for flat
             PARAMS[path[0]] = value;
        }
    };


    const createSliderRow = (label, uiKey, minVal, maxVal, stepVal, formatDigits = 2, ...paramPath) => {
        addParamToMapping(uiKey, ...paramPath);
        let currentValue = getParamValueFromPath(uiKey);
        if (currentValue === undefined) { console.error("Undefined PARAMS for UI (slider):", uiKey, paramPath); return; }

        let rowDiv = createDiv(label).parent(uiPanel).addClass('ui-row');
        let slider = createSlider(minVal, maxVal, currentValue, stepVal)
            .parent(rowDiv)
            .input(() => {
                let val = slider.value();
                setParamValueByPath(uiKey, val);
                select('#' + uiKey + 'ValueDisplay').html(nf(val, 0, formatDigits));
                if (!start) { initSimulation(); redraw(); }
            });
        createSpan(nf(currentValue, 0, formatDigits)).parent(rowDiv).id(uiKey + 'ValueDisplay');
        uiElements.sliders[uiKey] = slider;
        uiParamDisplays[uiKey] = select('#' + uiKey + 'ValueDisplay');
    };

    const createSelectRow = (label, uiKey, options, ...paramPath) => {
        addParamToMapping(uiKey, ...paramPath);
        let currentValue = getParamValueFromPath(uiKey);
        if (currentValue === undefined) { console.error("Undefined PARAMS for UI (select):", uiKey, paramPath); return; }

        let rowDiv = createDiv(label + ": ").parent(uiPanel).addClass('ui-row');
        let sel = createSelect().parent(rowDiv);
        options.forEach(opt => sel.option(typeof opt === 'string' ? opt : opt.label, typeof opt === 'string' ? opt : opt.value));
        sel.selected(currentValue);
        sel.changed(() => {
            setParamValueByPath(uiKey, sel.value());
            if (!start) { initSimulation(); redraw(); }
        });
        uiElements.selectors[uiKey] = sel;
    };
    const createSectionHeader = (text) => createP(`<b>${text}</b>`).parent(uiPanel);

    createSectionHeader("Simulation Controls:");
    let controlsDiv = createDiv().parent(uiPanel).id('controls-container');
    uiElements.buttons.startBtn = createButton("Start / Pause (Space)").parent(controlsDiv).mousePressed(toggleSimulation);
    uiElements.buttons.resetBtn = createButton("Reset Sim").parent(controlsDiv).mousePressed(() => { initSimulation(); redraw(); });
    uiElements.buttons.nutrientOverlayBtn = createButton("Show Nutrients").parent(controlsDiv).mousePressed(toggleNutrientOverlay);

    createSectionHeader("Environment:");
    createSliderRow("Nutrient Band Width:", 'W_band_slider', 5, 50, 1, 0, 'W_band');

    createSectionHeader("Cell & Strategy Settings:");
    createSelectRow("Initial Strategy Type:", 'default_strategy_type_sel', ["Responsive", "Constitutive", "BetHedging", "Memory"], 'default_strategy_type');
    createSliderRow("Initial Cluster Radius:", 'initial_cluster_radius_slider', 1, 10, 1, 0, 'initial_cluster_radius');
    createSliderRow("Initial Growth Rate G (lambda_G):", 'initial_inherent_growth_rate_G_slider', 0.1, PARAMS.max_possible_growth_G, 0.01, 2, 'initial_inherent_growth_rate_G');
    createSliderRow("Growth Rate L Factor (vs G):", 'lambda_L_factor_slider', 0.1, 1.5, 0.01, 2, 'lambda_L_factor');

    createSectionHeader("G->L Lag Trade-off & Modification:");
    createSliderRow("Base Lag (T_base):", 'trade_off_base_lag_slider', 0, 2.0, 0.05, 2, 'trade_off_base_lag');
    createSliderRow("Max Additional Lag (T_add_max):", 'trade_off_max_additional_lag_slider', 0, 5.0, 0.1, 1, 'trade_off_max_additional_lag');
    createSliderRow("Trade-off Exponent (k):", 'trade_off_exponent_k_slider', 0.5, 4.0, 0.1, 1, 'trade_off_exponent_k');
    createSliderRow("Parental Fitness Lag Mod Factor:", 'lag_mod_factor_slider', -2.0, 2.0, 0.05, 2, 'lag_mod_factor');

    createSectionHeader("Strategy: Constitutive");
    createSliderRow("Constitutive Growth G:", 'constitutive_growth_G_slider', 0.0, PARAMS.max_possible_growth_G, 0.01, 2, 'constitutive_growth_G');
    createSliderRow("Constitutive Growth L:", 'constitutive_growth_L_slider', 0.0, PARAMS.max_possible_growth_G, 0.01, 2, 'constitutive_growth_L');
    
    createSectionHeader("Strategy: BetHedging");
    createSliderRow("P(G parent -> L daughter):", 'bet_hedging_p_L_slider', 0, 1, 0.01, 2, 'bet_hedging_p_switch_to_L_if_G_parent');
    createSliderRow("P(L parent -> G daughter):", 'bet_hedging_p_G_slider', 0, 1, 0.01, 2, 'bet_hedging_p_switch_to_G_if_L_parent');

    createSectionHeader("Strategy: Memory");
    createSliderRow("Memory Window (time units):", 'memory_M_gen_memory_window_slider', 0.1, 20.0, 0.1, 1, 'memory_M_gen_memory_window');
    createSliderRow("Memory Lag Multiplier:", 'memory_T_lag_GL_multiplier_slider', 0.01, 1.0, 0.01, 2, 'memory_T_lag_GL_multiplier');


    createSectionHeader("Status:");
    let statusDiv = createDiv().parent(uiPanel).id('status-section');
    createSpan('Step: 0').parent(statusDiv).id('stepValueDisplay');
    createSpan('Time: 0.0').parent(statusDiv).id('timeValueDisplay');
    createSpan('Cell Count: 0').parent(statusDiv).id('cellCountValueDisplay');
    createSpan('Max Radius: 0.0').parent(statusDiv).id('maxRadiusValueDisplay');
    createSpan('Avg. Velocity: 0.00 u/t').parent(statusDiv).id('velocityValueDisplay');
    createSpan('Status: Initialized').parent(statusDiv).id('statusValueDisplay');

    updateUIParamDisplays(); // Populate initial values from PARAMS
}


function toggleSimulation() { start = !start; if (start) loop(); else noLoop(); updateDynamicUIText(); if (!start) redraw(); }
function toggleNutrientOverlay() { showNutrientOverlay = !showNutrientOverlay; uiElements.buttons.nutrientOverlayBtn.html(showNutrientOverlay ? 'Show Cells' : 'Show Nutrients'); if (!start) redraw(); }

function updateUIParamDisplays() { // Update displayed values in UI from PARAMS
    for (const uiKey in uiParamDisplays) {
        let valueToDisplay = getParamValueFromPath(uiKey); // Use the helper
        if (valueToDisplay === undefined) {
            // console.warn("Could not find param for UI display:", uiKey);
            continue;
        }

        let digits = 2; // Default digits for formatting
        const paramPath = paramMapping[uiKey];
        if (paramPath && paramPath.length > 0) {
            const actualKey = paramPath[paramPath.length -1]; // Get the final key in the path
            if (actualKey.includes('Rate') || actualKey.includes('Factor') || actualKey.includes('Multiplier') || actualKey.includes('Prob') || actualKey.includes('lag_mod')) digits = 2;
            else if (actualKey.includes('lag') || actualKey.includes('exponent') || actualKey.includes('Window') || actualKey.includes('max_additional')) digits = 1;
            else if (actualKey.includes('Radius') || actualKey.includes('Width') || actualKey.includes('band')) digits = 0;
        }
        
        if (uiParamDisplays[uiKey] && uiParamDisplays[uiKey].html() !== nf(valueToDisplay, 0, digits)) {
            uiParamDisplays[uiKey].html(nf(valueToDisplay, 0, digits));
        }
    }
}

function updateDynamicUIText() {
  select('#stepValueDisplay')?.html(`Step: ${simulationSteps}`);
  select('#timeValueDisplay')?.html(`Time: ${simulationTime.toFixed(2)}`);
  select('#cellCountValueDisplay')?.html(`Cell Count: ${liveCells.size}`);
  select('#maxRadiusValueDisplay')?.html(`Max Radius: ${maxColonyRadius.toFixed(2)}`);
  select('#velocityValueDisplay')?.html(`Avg. Velocity: ${avgRadialVelocity.toFixed(3)} u/t`);

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
// --- User Interaction ---
// =============================================================================
function mousePressed() {
  // Check if click is within the simulation buffer area on the canvas
  const bufferXOffset = (canvasWidth - bufferWidth) / 2;
  const bufferYOffset = (canvasHeight - bufferHeight) / 2;
  if (mouseX < bufferXOffset || mouseX > bufferXOffset + bufferWidth || 
      mouseY < bufferYOffset || mouseY > bufferYOffset + bufferHeight) {
      // Click was outside the simulation buffer area
      return;
  }
  // Convert mouse click on main canvas to buffer-relative coords
  let clickXInBuffer = mouseX - bufferXOffset;
  let clickYInBuffer = mouseY - bufferYOffset;
  
  // Convert buffer-relative coords to grid-relative coords for pixelToGridCoords
  // pixelToGridCoords expects coords relative to grid (0,0) pixel space which is gridOriginX/Y from buffer's top-left
  const gridCoords = pixelToGridCoords(clickXInBuffer, clickYInBuffer); 
  const coordStr = `${gridCoords.x},${gridCoords.y}`;

  if (liveCells.has(coordStr)) {
    console.log(`Cell at (${gridCoords.x}, ${gridCoords.y}):`, liveCells.get(coordStr));
  } else {
    console.log(`Empty slot at (${gridCoords.x}, ${gridCoords.y}).`);
  }
  if (nutrientGrid.has(coordStr)) {
    console.log(`Nutrient: ${nutrientGrid.get(coordStr)}`);
  } else {
    console.log(`Nutrient: Out of bounds.`);
  }
}

function keyPressed() { if (key === ' ') { toggleSimulation(); return false; } }
// =============================================================================
