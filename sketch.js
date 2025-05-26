// --- Core Microbial Strategy Simulation (v32 - Square Lattice, Concentric Bands, Strategy Types) ---
// Model based on "2D Radial Microbial Range Expansion in Fluctuating Environments" description
// Key Features:
// - Square Lattice Environment
// - Concentric Alternating Nutrient Bands (Glucose/Galactose)
// - Cell Agents with Strategies: Constitutive, Responsive, BetHedging, Memory
// - Growth Rate vs. Lag Time Trade-off
// - DK-Inspired Lag Modification for New Cells
// - Discrete Time Steps (dt)
// ---

// --- Simulation Settings ---
let worldSize = 120; // Half-width of the square world (e.g., -worldSize to +worldSize)
let cellSize = 5;    // Pixel size of each cell
let dt = 0.1;        // Simulation time step
const epsilon = 1e-9;
let drawInterval = 1;

// --- MODEL PARAMETERS (Defaults, adjustable via UI) ---
let PARAMS = {
  // --- Environment ---
  W_band: 15, // Width of each nutrient band (in lattice units)

  // --- Cell General ---
  initial_cluster_radius: 3, // Radius of initial G-specialist cell cluster
  default_strategy_type: "Responsive", // Strategy for initial cells
  lambda_L_factor: 0.7, // Growth rate on Galactose = lambda_L_factor * inherent_growth_rate_G
  max_possible_growth_G: 1.0, // For normalization in trade-offs/lag calcs

  // --- Trade-off: Growth Rate G vs. Lag Time G->L ---
  trade_off_base_lag: 0.5,
  trade_off_max_additional_lag: 3.0, // Max T_lag_GL will be base_lag + this
  trade_off_exponent_k: 2.0, // Curvature of trade-off

  // --- DK-Inspired Lag Modification for New Cells ---
  // Factor for inherent_T_lag_GL_daughter * (1 + lag_mod_factor * (avg_fitness_G_parents / max_possible_growth_G))
  lag_mod_factor: 0.5, // Positive: better G-parents increase daughter's G->L lag. Negative would decrease.

  // --- Strategy: Constitutive ---
  constitutive_growth_G: 0.6, // Growth rate on G for constitutive
  constitutive_growth_L: 0.5, // Growth rate on L for constitutive

  // --- Strategy: BetHedging ---
  bet_hedging_p_switch_to_L_if_G_parent: 0.1, // Prob G-parent makes L-daughter
  bet_hedging_p_switch_to_G_if_L_parent: 0.1, // Prob L-parent makes G-daughter

  // --- Strategy: Memory ---
  memory_M_gen_memory_window: 5.0,  // Time units (e.g. 5 * dt) cell remembers being on Gal
  memory_T_lag_GL_multiplier: 0.2, // Lag is X% of inherent if within memory window

  // --- (Future: parameters for evolving strategies, if added) ---
  // inherent_growth_rate_G will be a property of each cell lineage
  // For now, we can use a default or allow UI to set for initial population.
  initial_inherent_growth_rate_G: 0.8,
};

// --- Simulation State ---
let nutrientGrid = new Map(); // Stores "Glucose" or "Galactose" for each (x,y) string key "x,y"
let liveCells = new Map();    // Stores cell objects, key is "x,y" string
let nextCellID = 0;
let simulationTime = 0;
let simulationSteps = 0;
let maxColonyRadius = 0;      // Max Euclidean distance of a cell from origin
let lastStepBirths = 1;       // To detect stasis
let start = false;
let showNutrientOverlay = false;

// --- Graphics & Grid ---
let nutrientBuffer, cellBuffer; // cellBuffer for drawing cells efficiently
let finalStateBuffer; let finalStateRendered = false;
let canvasWidth, canvasHeight, bufferWidth, bufferHeight;
let gridOriginX, gridOriginY; // Pixel coords of (0,0) in the grid

// Moore neighborhood (8 neighbors) for finding adjacent cells/slots
const moore_deltas = [
  { dx: -1, dy: -1 }, { dx: 0, dy: -1 }, { dx: 1, dy: -1 },
  { dx: -1, dy: 0 },  /* {dx:0, dy:0} */ { dx: 1, dy: 0 },
  { dx: -1, dy: 1 },  { dx: 0, dy: 1 },  { dx: 1, dy: 1 },
];
// Von Neumann neighborhood (4 neighbors), if needed for specific interactions later
// const von_neumann_deltas = [ {dx:0, dy:1}, {dx:0, dy:-1}, {dx:1, dy:0}, {dx:-1, dy:0} ];

// --- Color Definitions ---
let colorGlucoseBand, colorGalactoseBand;
let colorGSpecialist, colorLSpecialist, colorSwitchingGL;
let strategyColors = {}; // For potential visual distinction of strategies

// --- UI Elements ---
let uiElements = { controlGroups: {} };
let uiParamDisplays = {};


// =============================================================================
// --- Math & Grid Helpers (NEW - Square Lattice) ---
// =============================================================================
function euclideanDistance(x1, y1, x2, y2) {
  return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
}

function pixelToGridCoords(pixelX, pixelY) {
  // Convert screen pixel coordinates to integer grid coordinates
  // Assumes gridOriginX, gridOriginY are the pixel coordinates of lattice (0,0)
  let worldX = pixelX - gridOriginX;
  let worldY = pixelY - gridOriginY;
  let gridX = Math.floor(worldX / cellSize);
  let gridY = Math.floor(worldY / cellSize);
  return { x: gridX, y: gridY };
}

function gridToPixelCoords(gridX, gridY) {
  // Convert integer grid coordinates to top-left pixel of the cell
  // Relative to gridOriginX, gridOriginY
  return { x: gridX * cellSize, y: gridY * cellSize };
}

// =============================================================================
// --- Drawing Helpers (NEW - Square Lattice) ---
// =============================================================================
function calculateGridDimensions() {
  bufferWidth = (worldSize * 2 + 1) * cellSize;
  bufferHeight = (worldSize * 2 + 1) * cellSize;
  canvasWidth = bufferWidth + 40; // some padding
  canvasHeight = bufferHeight + 40; // some padding
  gridOriginX = canvasWidth / 2;
  gridOriginY = canvasHeight / 2;
}

function drawSquareCell(buffer, gridX, gridY, fillColor, strokeColor = null, weight = 0.5) {
  const p = gridToPixelCoords(gridX, gridY); // Gets top-left relative to origin
  buffer.fill(fillColor);
  if (strokeColor) {
    buffer.stroke(strokeColor);
    buffer.strokeWeight(weight);
  } else {
    buffer.noStroke();
  }
  // Draw relative to the buffer's origin, which is top-left of buffer.
  // The gridToPixelCoords gives coords relative to (0,0) of the grid.
  // When drawing on a buffer translated by gridOrigin (which matches buffer center),
  // these coords are correct.
  buffer.rect(p.x, p.y, cellSize, cellSize);
}


// =============================================================================
// --- P5.js Core Functions ---
// =============================================================================
function setup() {
  calculateGridDimensions();
  nutrientBuffer = createGraphics(bufferWidth, bufferHeight);
  cellBuffer = createGraphics(bufferWidth, bufferHeight);
  finalStateBuffer = createGraphics(bufferWidth, bufferHeight);
  finalStateBuffer.pixelDensity(1); // For saving final state

  console.log(`Grid/Buffer Size: ${bufferWidth}x${bufferHeight}. Cell Size: ${cellSize}`);

  // Nutrient Environment Colors
  colorGlucoseBand = color(255, 230, 230, 180); // Light reddish/pink for Glucose band
  colorGalactoseBand = color(230, 230, 255, 180); // Light bluish for Galactose band

  // Cell Phenotype Colors
  colorGSpecialist = color(0, 200, 0);     // Green
  colorLSpecialist = color(0, 0, 200);     // Blue
  colorSwitchingGL = color(255, 165, 0);   // Orange

  // (Optional) Strategy Colors - could be used to tint phenotype color or draw an outline
  strategyColors["Constitutive"] = color(128, 128, 128, 100); // Grey
  strategyColors["Responsive"] = color(200, 200, 0, 100);   // Yellowish
  strategyColors["BetHedging"] = color(200, 0, 200, 100);   // Magentaish
  strategyColors["Memory"] = color(0, 200, 200, 100);     // Cyanish

  let totalCanvasHeight = Math.max(canvasHeight, 750);
  let mainCanvas = createCanvas(canvasWidth, totalCanvasHeight);
  let canvasContainer = select('#canvas-container');
  if (canvasContainer) mainCanvas.parent('canvas-container');
  else console.error("#canvas-container not found.");

  createUI(); // Will be significantly different
  initSimulation();

  background(240);
  push();
  translate(gridOriginX, gridOriginY); // Translate once for drawing buffers
  image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2);
  image(cellBuffer, -bufferWidth / 2, -bufferHeight / 2);
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
    translate(gridOriginX, gridOriginY); // Center the drawing space

    if (showNutrientOverlay) {
      image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2);
    } else {
      // Draw nutrient buffer first if not in overlay mode, then cells
      image(nutrientBuffer, -bufferWidth / 2, -bufferHeight / 2); // Always show nutrients underneath
      if (!start && finalStateRendered) {
        image(finalStateBuffer, -bufferWidth / 2, -bufferHeight / 2);
      } else {
        drawLiveCellsOnBuffer(cellBuffer); // Redraw cells to their buffer
        image(cellBuffer, -bufferWidth / 2, -bufferHeight / 2);
      }
    }
    pop();
  }
  updateDynamicUIText();
}

// =============================================================================
// --- Nutrient Environment Logic (NEW - Concentric Bands) ---
// =============================================================================
function getNutrientAt(x, y) { // Grid coordinates
  const dist = euclideanDistance(x, y, 0, 0);
  const band_index = Math.floor(dist / PARAMS.W_band);
  if (band_index % 2 === 0) {
    return "Glucose"; // Even bands (0, 2, 4...) are Glucose
  } else {
    return "Galactose"; // Odd bands (1, 3, 5...) are Galactose
  }
}

function populateNutrientGridAndDrawBuffer() {
  nutrientGrid.clear();
  nutrientBuffer.push();
  nutrientBuffer.translate(bufferWidth / 2, bufferHeight / 2); // Center drawing in buffer
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
  if (start) redraw();
}


// =============================================================================
// --- Cell Agent Logic (NEW) ---
// =============================================================================

// Key function: Trade-off between growth rate on G and G->L lag time
function trade_off_lag_vs_growth(inherent_growth_rate_G) {
  const { trade_off_base_lag, trade_off_max_additional_lag, trade_off_exponent_k, max_possible_growth_G } = PARAMS;
  if (inherent_growth_rate_G <= epsilon) return trade_off_base_lag + trade_off_max_additional_lag; // Max lag if no G growth
  
  // Ensure growth rate doesn't exceed max_possible_growth_G for this calculation
  const normalized_gr_G = Math.min(inherent_growth_rate_G, max_possible_growth_G) / max_possible_growth_G;
  
  let lag = trade_off_base_lag + Math.pow(normalized_gr_G, trade_off_exponent_k) * trade_off_max_additional_lag;
  return Math.max(0, lag); // Ensure lag is not negative
}

// Key function: DK-Inspired lag modification for new cells
function calculate_modified_lag(inherent_T_lag_GL_daughter, avg_fitness_G_parents_of_slot) {
  const { lag_mod_factor, max_possible_growth_G } = PARAMS;
  if (avg_fitness_G_parents_of_slot < epsilon) return inherent_T_lag_GL_daughter; // No G-fit parents, no modification

  // Normalize avg_fitness for this calculation
  const normalized_avg_fitness = Math.min(avg_fitness_G_parents_of_slot, max_possible_growth_G) / max_possible_growth_G;
  
  let modifiedLag = inherent_T_lag_GL_daughter * (1 + lag_mod_factor * normalized_avg_fitness);
  return Math.max(0, modifiedLag); // Ensure lag is not negative
}

function createCell(x, y, strategy_type, inherent_growth_rate_G, parentCell = null) {
  nextCellID++;
  const inherent_T_lag_GL = trade_off_lag_vs_growth(inherent_growth_rate_G);
  
  let newCell = {
    id: nextCellID,
    x: x, y: y,
    phenotype: "G_specialist", // Default, will be adjusted
    strategy_type: strategy_type,
    inherent_growth_rate_G: inherent_growth_rate_G,
    inherent_T_lag_GL: inherent_T_lag_GL,
    current_growth_rate: 0.0,
    remaining_lag_time: 0.0,
    memory_state: 0.0, // Time since last on Galactose for "Memory" strategy
    // Could add parent_id if lineage tracking is needed
  };

  if (parentCell) { // Inherit memory state if applicable
      newCell.memory_state = parentCell.memory_state; // Or specific logic for memory inheritance
  }
  
  // Initial phenotype and growth rate will be set after creation based on local nutrient
  // and potentially parent/bet-hedging logic.
  return newCell;
}


// =============================================================================
// --- Initialization & Reset ---
// =============================================================================
function initSimulation() {
  start = false; loop();
  lastStepBirths = 1; finalStateRendered = false; showNutrientOverlay = false;
  if (uiElements.buttons?.nutrientOverlayBtn) uiElements.buttons.nutrientOverlayBtn.html('Show Nutrients');
  if (finalStateBuffer) finalStateBuffer.clear();
  cellBuffer.clear(); // Clear cell drawing buffer

  simulationTime = 0;
  simulationSteps = 0;
  maxColonyRadius = 0;
  nextCellID = 0;

  liveCells.clear();
  populateNutrientGridAndDrawBuffer(); // Sets up nutrientGrid and nutrientBuffer

  // Initialize a small, dense circular cluster of G_specialist cells
  const { initial_cluster_radius, default_strategy_type, initial_inherent_growth_rate_G } = PARAMS;
  let plantedCount = 0;
  for (let gx = -initial_cluster_radius; gx <= initial_cluster_radius; gx++) {
    for (let gy = -initial_cluster_radius; gy <= initial_cluster_radius; gy++) {
      if (euclideanDistance(gx, gy, 0, 0) <= initial_cluster_radius) {
        if (nutrientGrid.has(`${gx},${gy}`)) { // Ensure it's within world bounds
          const coordStr = `${gx},${gy}`;
          if (!liveCells.has(coordStr)) {
            let cell = createCell(gx, gy, default_strategy_type, initial_inherent_growth_rate_G);
            
            // Initial state update for these first cells
            const local_nutrient = nutrientGrid.get(coordStr);
            if (local_nutrient === "Glucose") {
                cell.phenotype = "G_specialist";
                cell.current_growth_rate = cell.inherent_growth_rate_G;
                 // For Memory strategy: if starting on Glucose, they haven't seen Gal yet
                if (cell.strategy_type === "Memory") cell.memory_state = PARAMS.memory_M_gen_memory_window + 1; // Effectively no recent memory
            } else { // Started on Galactose - this should be rare if W_band is large enough for initial cluster
                cell.phenotype = "Switching_GL";
                cell.remaining_lag_time = cell.inherent_T_lag_GL; // No parental modification for initial seed
                cell.current_growth_rate = 0.0;
                 if (cell.strategy_type === "Memory") cell.memory_state = 0; // Starting on Gal, so fresh memory
            }
            liveCells.set(coordStr, cell);
            maxColonyRadius = Math.max(maxColonyRadius, euclideanDistance(gx, gy, 0, 0));
            plantedCount++;
          }
        }
      }
    }
  }
  console.log(`Initialized ${plantedCount} cells with strategy ${default_strategy_type}. Max Radius: ${maxColonyRadius.toFixed(2)}`);
  
  drawLiveCellsOnBuffer(cellBuffer); // Initial draw
  noLoop(); // Start paused
  updateDynamicUIText();
}

// =============================================================================
// --- Simulation Step Logic (NEW) ---
// =============================================================================
function stepSimulation() {
  if (liveCells.size === 0 || (lastStepBirths === 0 && simulationSteps > 0)) {
    if (start) {
      let reason = liveCells.size === 0 ? "Extinction" : "Stasis";
      console.log(`${reason} @ step ${simulationSteps}, t=${simulationTime.toFixed(2)}. Stopping.`);
      if (liveCells.size > 0 && !finalStateRendered) {
        drawLiveCellsOnBuffer(finalStateBuffer); // Use cellBuffer logic for final state
        finalStateRendered = true;
      }
      start = false; noLoop(); updateDynamicUIText();
    }
    return;
  }

  simulationSteps++;
  simulationTime += dt;
  let currentBirthsThisStep = 0;

  // --- 1. Cell State Update (for existing cells) ---
  let cellsToUpdate = Array.from(liveCells.values()); // Create a copy to iterate over, as liveCells might change if we had death
  
  for (const cell of cellsToUpdate) {
    const local_nutrient = nutrientGrid.get(`${cell.x},${cell.y}`);
    let previous_phenotype = cell.phenotype;

    // A. Lag countdown
    if (cell.phenotype === "Switching_GL") {
      cell.remaining_lag_time -= dt;
      if (cell.remaining_lag_time <= epsilon) {
        cell.remaining_lag_time = 0;
        cell.phenotype = "L_specialist";
        if (cell.strategy_type === "Memory") cell.memory_state = 0; // Just switched to L, reset memory timer
      }
    }
    // B. Strategy-specific phenotype switching (if not lagging)
    else {
        // --- Responsive Strategy ---
        if (cell.strategy_type === "Responsive") {
            if (cell.phenotype === "G_specialist" && local_nutrient === "Galactose") {
                cell.phenotype = "Switching_GL";
                cell.remaining_lag_time = cell.inherent_T_lag_GL;
                // Memory modification for Responsive strategy when it initiates lag
                if (PARAMS.memory_M_gen_memory_window > 0 && cell.memory_state < PARAMS.memory_M_gen_memory_window) {
                    cell.remaining_lag_time *= PARAMS.memory_T_lag_GL_multiplier;
                }
            } else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") {
                cell.phenotype = "G_specialist"; // No L->G lag
            }
        }
        // --- Memory Strategy ---
        else if (cell.strategy_type === "Memory") {
            if (cell.phenotype === "G_specialist" && local_nutrient === "Galactose") {
                cell.phenotype = "Switching_GL";
                // Apply memory: if recently on Gal (low memory_state), reduce lag
                if (cell.memory_state < PARAMS.memory_M_gen_memory_window) {
                    cell.remaining_lag_time = cell.inherent_T_lag_GL * PARAMS.memory_T_lag_GL_multiplier;
                } else {
                    cell.remaining_lag_time = cell.inherent_T_lag_GL;
                }
            } else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") {
                cell.phenotype = "G_specialist";
            }
        }
        // --- Constitutive Strategy --- (No phenotype switching, growth rate handles it)
        // --- BetHedging Strategy --- (Switching handled at birth, existing cells are responsive-like)
        else if (cell.strategy_type === "BetHedging") { // Behaves like responsive once established
             if (cell.phenotype === "G_specialist" && local_nutrient === "Galactose") {
                cell.phenotype = "Switching_GL";
                cell.remaining_lag_time = cell.inherent_T_lag_GL;
            } else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") {
                cell.phenotype = "G_specialist";
            }
        }
    }

    // C. Update current_growth_rate
    cell.current_growth_rate = 0.0; // Default, esp. if lagging
    if (cell.phenotype === "Switching_GL") {
      // Growth rate is 0 while switching
    } else if (cell.strategy_type === "Constitutive") {
        cell.current_growth_rate = (local_nutrient === "Glucose") ? PARAMS.constitutive_growth_G : PARAMS.constitutive_growth_L;
    } else { // For Responsive, Memory, BetHedging (once phenotype set)
        if (cell.phenotype === "G_specialist") {
            if (local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G;
            // else (G_specialist on Galactose, but not yet switched to Switching_GL, e.g. first frame) growth is 0
        } else if (cell.phenotype === "L_specialist") {
            if (local_nutrient === "Galactose") cell.current_growth_rate = PARAMS.lambda_L_factor * cell.inherent_growth_rate_G;
            else if (local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G; // L on G, grows fine
        }
    }
    
    // D. Update memory_state (for Memory & Responsive if using memory effect)
    if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
        if (cell.phenotype === "L_specialist" && local_nutrient === "Galactose") {
            cell.memory_state = 0; // Reset timer: on Galactose
        } else if (cell.phenotype === "Switching_GL" && local_nutrient === "Galactose") {
             // If switching *to* L, consider memory active. Or reset if switch completes.
             // This is handled when switch to L_specialist completes.
        }
        else {
            cell.memory_state += dt; // Increment time since last definitively on Gal
        }
    }
  } // End of cell state update loop

  // --- 2. New Cell Generation (at Frontier) ---
  let frontierSlots = new Map(); // Key: "x,y" of empty slot, Value: {parentalNeighbors: [cellObj, ...]}

  for (const [coordStr, cell] of liveCells) {
    if (cell.current_growth_rate < epsilon && cell.phenotype !== "Switching_GL") continue; // Skip if not growing and not just lagging

    for (const delta of moore_deltas) {
      const nx = cell.x + delta.dx;
      const ny = cell.y + delta.dy;
      const nCoordStr = `${nx},${ny}`;

      if (nx >= -worldSize && nx <= worldSize && ny >= -worldSize && ny <= worldSize) { // Check bounds
        if (!liveCells.has(nCoordStr)) { // If slot is empty
          if (!frontierSlots.has(nCoordStr)) {
            frontierSlots.set(nCoordStr, { parentalNeighbors: [] });
          }
          // Add current cell `cell` as a potential parent for this empty slot `nCoordStr`
          // We will collect all potential parents for each slot first
        }
      }
    }
  }
  
  // Now, populate parentalNeighbors for each identified frontier slot
  for (const [slotCoordStr, slotData] of frontierSlots) {
      const [sx, sy] = slotCoordStr.split(',').map(Number);
      for (const delta of moore_deltas) { // Check neighbors of the slot
          const px = sx + delta.dx;
          const py = sy + delta.dy;
          const pCoordStr = `${px},${py}`;
          if (liveCells.has(pCoordStr)) {
              slotData.parentalNeighbors.push(liveCells.get(pCoordStr));
          }
      }
  }


  let newCellsToAdd = []; // Store {coordStr: "x,y", cell: newCellObject}

  for (const [slotCoordStr, slotData] of frontierSlots) {
    if (slotData.parentalNeighbors.length === 0) continue;

    // Stochastic Colonization: Choose parent
    let chosenParent = null;
    let totalGrowthRate = 0;
    slotData.parentalNeighbors.forEach(p => totalGrowthRate += p.current_growth_rate);

    if (totalGrowthRate > epsilon) {
      let randVal = random() * totalGrowthRate;
      let cumulativeGrowth = 0;
      for (const p of slotData.parentalNeighbors) {
        cumulativeGrowth += p.current_growth_rate;
        if (randVal <= cumulativeGrowth) {
          chosenParent = p;
          break;
        }
      }
    } else { // All potential parents have 0 growth rate (e.g., all lagging)
      chosenParent = random(slotData.parentalNeighbors);
    }

    if (!chosenParent) continue; // Should not happen if parentalNeighbors > 0

    const [nx, ny] = slotCoordStr.split(',').map(Number);
    let daughterCell = createCell(nx, ny, chosenParent.strategy_type, chosenParent.inherent_growth_rate_G, chosenParent);
    
    // Daughter's Initial Phenotype and Lag
    const local_nutrient_at_slot = nutrientGrid.get(slotCoordStr);
    daughterCell.phenotype = chosenParent.phenotype; // Initial inheritance

    // BetHedging specific logic for initial phenotype
    if (daughterCell.strategy_type === "BetHedging") {
        if (chosenParent.phenotype === "G_specialist" && random() < PARAMS.bet_hedging_p_switch_to_L_if_G_parent) {
            daughterCell.phenotype = "L_specialist"; // Bet-hedged to L
        } else if (chosenParent.phenotype === "L_specialist" && random() < PARAMS.bet_hedging_p_switch_to_G_if_L_parent) {
            daughterCell.phenotype = "G_specialist"; // Bet-hedged to G
        }
    }
    
    // Lag Modification for G->L switch
    if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Galactose") {
        daughterCell.phenotype = "Switching_GL";
        let avg_fitness_G_parents = 0;
        let g_parent_count = 0;
        slotData.parentalNeighbors.forEach(p => {
            // Consider 'fitness' for DK-lag to be inherent_growth_rate_G
            // Or could be current_growth_rate if on G. Here using inherent_growth_rate_G as per description hint.
            if (p.phenotype === "G_specialist" || p.strategy_type === "Constitutive") { // Constitutive always has some G potential
                 avg_fitness_G_parents += p.inherent_growth_rate_G; // Or p.current_growth_rate if they are on G
                 g_parent_count++;
            }
        });
        if (g_parent_count > 0) avg_fitness_G_parents /= g_parent_count;
        else avg_fitness_G_parents = 0; // Or some default low fitness

        daughterCell.remaining_lag_time = calculate_modified_lag(daughterCell.inherent_T_lag_GL, avg_fitness_G_parents);

        // Memory consideration for daughter IF its strategy is Memory or Responsive(with memory)
        if (daughterCell.strategy_type === "Memory" || (daughterCell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0) ) {
            if (chosenParent.memory_state < PARAMS.memory_M_gen_memory_window) { // Parent had recent Gal memory
                 daughterCell.remaining_lag_time = Math.min(daughterCell.remaining_lag_time, daughterCell.inherent_T_lag_GL * PARAMS.memory_T_lag_GL_multiplier);
            }
        }

    } else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Glucose") {
        daughterCell.phenotype = "G_specialist"; // L->G switch, no lag
        daughterCell.remaining_lag_time = 0;
    } else { // No switch needed or L->G (which is instant)
        daughterCell.remaining_lag_time = 0;
    }

    // Set daughter's initial growth rate based on its final initial phenotype and local nutrient
    daughterCell.current_growth_rate = 0.0; // Default
    if (daughterCell.phenotype === "Switching_GL") {
        // Growth is 0
    } else if (daughterCell.strategy_type === "Constitutive") {
        daughterCell.current_growth_rate = (local_nutrient_at_slot === "Glucose") ? PARAMS.constitutive_growth_G : PARAMS.constitutive_growth_L;
    } else {
        if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Glucose") {
            daughterCell.current_growth_rate = daughterCell.inherent_growth_rate_G;
        } else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Galactose") {
            daughterCell.current_growth_rate = PARAMS.lambda_L_factor * daughterCell.inherent_growth_rate_G;
        } else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Glucose") { // Already switched to G_specialist if needed
            daughterCell.current_growth_rate = daughterCell.inherent_growth_rate_G;
        }
    }
    
    // Update daughter's memory state based on initial placement
    if (daughterCell.strategy_type === "Memory" || (daughterCell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
        if ((daughterCell.phenotype === "L_specialist" || daughterCell.phenotype === "Switching_GL") && local_nutrient_at_slot === "Galactose") {
            daughterCell.memory_state = 0; // Born onto Gal, fresh memory
        } else { // Born onto Glu, or parent was on Glu
            daughterCell.memory_state = chosenParent.memory_state + dt; // Inherit and advance, or just advance if parent memory not applicable
        }
    }

    newCellsToAdd.push({ coordStr: slotCoordStr, cell: daughterCell });
    currentBirthsThisStep++;
    maxColonyRadius = Math.max(maxColonyRadius, euclideanDistance(nx, ny, 0, 0));
  } // End of frontier slot loop

  // Add all new cells to liveCells
  newCellsToAdd.forEach(item => {
    liveCells.set(item.coordStr, item.cell);
  });

  lastStepBirths = currentBirthsThisStep;
  if (liveCells.size === 0) lastStepBirths = 0; // Ensure stasis if extinct
}


// =============================================================================
// --- Drawing Functions (Cell drawing NEW) ---
// =============================================================================
function drawLiveCellsOnBuffer(targetBuffer) {
  if (!targetBuffer) return;
  targetBuffer.clear(); // Clear previous frame on this buffer
  targetBuffer.push();
  targetBuffer.translate(targetBuffer.width / 2, targetBuffer.height / 2); // Center drawing

  for (const [coordStr, cell] of liveCells.entries()) {
    let cellCol;
    if (cell.phenotype === "Switching_GL") {
      cellCol = colorSwitchingGL;
    } else if (cell.phenotype === "G_specialist") {
      cellCol = colorGSpecialist;
    } else { // L_specialist
      cellCol = colorLSpecialist;
    }
    // Optional: Add a slight tint or outline based on strategy_type if desired
    // e.g., lerpColor(cellCol, strategyColors[cell.strategy_type], 0.3)
    drawSquareCell(targetBuffer, cell.x, cell.y, cellCol, color(50), 0.2);
  }
  targetBuffer.pop();
}


// =============================================================================
// --- UI Setup & Updates (Needs MAJOR rework) ---
// =============================================================================
function createUI() {
    const uiPanel = select('#ui-panel'); uiPanel.html('');
    uiElements = { controlGroups: {}, sliders: {}, selectors: {}, buttons: {} };
    uiParamDisplays = {};

    const createSliderRow = (label, paramKey, minVal, maxVal, stepVal, formatDigits = 2, subKey = null, requiresReset = true) => {
        let path = subKey ? PARAMS[paramKey][subKey] : PARAMS[paramKey];
        if (path === undefined && subKey) path = PARAMS[subKey]; // Direct key if not nested
        else if (path === undefined) { console.error("Undefined PARAMS path for UI:", paramKey, subKey); return; }


        let rowDiv = createDiv(label).parent(uiPanel).addClass('ui-row');
        let slider = createSlider(minVal, maxVal, path, stepVal)
            .parent(rowDiv)
            .input(() => {
                let val = slider.value();
                if (subKey && PARAMS[paramKey] && typeof PARAMS[paramKey] === 'object') PARAMS[paramKey][subKey] = val;
                else PARAMS[paramKey] = val; // Fallback for direct key or if subKey logic was simpler before

                select('#' + (subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay').html(nf(val, 0, formatDigits));
                if (requiresReset && !start) {
                    console.log(`Param ${paramKey}${subKey ? '.'+subKey : ''} changed, re-init simulation...`);
                    initSimulation(); redraw();
                }
            });
        createSpan(nf(path, 0, formatDigits)).parent(rowDiv).id((subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay');
        uiElements.sliders[(subKey ? paramKey + '_' + subKey : paramKey)] = slider;
        uiParamDisplays[(subKey ? paramKey + '_' + subKey : paramKey)] = select('#' + (subKey ? paramKey + '_' + subKey : paramKey) + 'ValueDisplay');
        return rowDiv;
    };
    const createSectionHeader = (text) => createP(`<b>${text}</b>`).parent(uiPanel);
    const createSelectRow = (label, paramKey, options, requiresReset = true) => {
        let rowDiv = createDiv(label + ": ").parent(uiPanel).addClass('ui-row');
        let sel = createSelect().parent(rowDiv);
        options.forEach(opt => sel.option(typeof opt === 'string' ? opt : opt.label, typeof opt === 'string' ? opt : opt.value));
        sel.selected(PARAMS[paramKey]);
        sel.changed(() => {
            PARAMS[paramKey] = sel.value();
            console.log(`${paramKey} changed to ${sel.value()}`);
            if (requiresReset && !start) {
                initSimulation(); redraw();
            }
        });
        uiElements.selectors[paramKey] = sel; return rowDiv;
    };

    createSectionHeader("Simulation Controls:");
    let controlsDiv = createDiv().parent(uiPanel).id('controls-container');
    uiElements.buttons.startBtn = createButton("Start / Pause (Space)").parent(controlsDiv).mousePressed(toggleSimulation);
    uiElements.buttons.resetBtn = createButton("Reset Sim").parent(controlsDiv).mousePressed(() => { initSimulation(); redraw(); });
    uiElements.buttons.nutrientOverlayBtn = createButton("Show Nutrients").parent(controlsDiv).mousePressed(toggleNutrientOverlay);

    createSectionHeader("Environment:");
    createSliderRow("Nutrient Band Width (W_band):", 'W_band', 5, 50, 1, 0);

    createSectionHeader("Cell & Strategy Settings:");
    createSelectRow("Initial Strategy Type:", 'default_strategy_type', ["Responsive", "Constitutive", "BetHedging", "Memory"]);
    createSliderRow("Initial Cluster Radius:", 'initial_cluster_radius', 1, 10, 1, 0);
    createSliderRow("Initial Growth Rate G (lambda_G):", 'initial_inherent_growth_rate_G', 0.1, PARAMS.max_possible_growth_G, 0.01, 2);
    createSliderRow("Growth Rate L Factor (vs G):", 'lambda_L_factor', 0.1, 1.5, 0.01, 2);

    createSectionHeader("G->L Lag Trade-off & Modification:");
    createSliderRow("Base Lag (T_base):", 'trade_off_base_lag', 0, 2.0, 0.05, 2);
    createSliderRow("Max Additional Lag (T_add_max):", 'trade_off_max_additional_lag', 0, 5.0, 0.1, 1);
    createSliderRow("Trade-off Exponent (k):", 'trade_off_exponent_k', 0.5, 4.0, 0.1, 1);
    createSliderRow("Parental Fitness Lag Mod Factor:", 'lag_mod_factor', -1.0, 1.0, 0.05, 2);

    createSectionHeader("Strategy: Constitutive");
    createSliderRow("Constitutive Growth G:", 'constitutive_growth_G', 0.0, PARAMS.max_possible_growth_G, 0.01, 2);
    createSliderRow("Constitutive Growth L:", 'constitutive_growth_L', 0.0, PARAMS.max_possible_growth_G, 0.01, 2);
    
    createSectionHeader("Strategy: BetHedging");
    createSliderRow("P(G parent -> L daughter):", 'bet_hedging_p_switch_to_L_if_G_parent', 0, 1, 0.01, 2);
    createSliderRow("P(L parent -> G daughter):", 'bet_hedging_p_switch_to_G_if_L_parent', 0, 1, 0.01, 2);

    createSectionHeader("Strategy: Memory");
    createSliderRow("Memory Window (time units):", 'memory_M_gen_memory_window', 0.1, 20.0, 0.1, 1);
    createSliderRow("Memory Lag Multiplier:", 'memory_T_lag_GL_multiplier', 0.01, 1.0, 0.01, 2);


    createSectionHeader("Status:");
    let statusDiv = createDiv().parent(uiPanel).id('status-section');
    createSpan('Step: 0').parent(statusDiv).id('stepValueDisplay');
    createSpan('Time: 0.0').parent(statusDiv).id('timeValueDisplay');
    createSpan('Cell Count: 0').parent(statusDiv).id('cellCountValueDisplay');
    createSpan('Max Radius: 0.0').parent(statusDiv).id('maxRadiusValueDisplay');
    // TODO: Add v_rad display
    createSpan('Status: Initialized').parent(statusDiv).id('statusValueDisplay');

    updateUIParamDisplays();
}


function toggleSimulation() { start = !start; if (start) loop(); else noLoop(); updateDynamicUIText(); if (!start) redraw(); }
function toggleNutrientOverlay() { showNutrientOverlay = !showNutrientOverlay; uiElements.buttons.nutrientOverlayBtn.html(showNutrientOverlay ? 'Show Cells' : 'Show Nutrients'); if (!start) redraw(); }

function updateUIParamDisplays() {
    for (const key in uiParamDisplays) {
        let valueToDisplay, pKey = key, subKey = null, digits = 2; // Default 2 digits
        if (key.includes('_')) { // Simple check for nested, might need better logic if deeply nested keys
            const parts = key.split('_'); // This splitting is naive for complex keys
            pKey = parts[0];
            // Attempt to reconstruct subKey, assuming first part is main key
            if (PARAMS[pKey] && typeof PARAMS[pKey] === 'object') {
                 subKey = key.substring(pKey.length + 1);
                 if (PARAMS[pKey].hasOwnProperty(subKey)) valueToDisplay = PARAMS[pKey][subKey];
                 else valueToDisplay = PARAMS[key]; // Fallback: key is direct PARAMS property
            } else {
                valueToDisplay = PARAMS[key]; // Key is direct PARAMS property
            }
        } else if (PARAMS.hasOwnProperty(key)) {
            valueToDisplay = PARAMS[key];
        } else continue;

        if (valueToDisplay === undefined) continue;

        // Infer digits (very simplified for this version)
        if (key.includes('Rate') || key.includes('Factor') || key.includes('Multiplier') || key.includes('Prob')) digits = 2;
        else if (key.includes('lag') || key.includes('exponent') || key.includes('Window')) digits = 1;
        else if (key.includes('Radius') || key.includes('Width')) digits = 0;
        
        if (uiParamDisplays[key] && uiParamDisplays[key].html() !== nf(valueToDisplay, 0, digits)) {
            uiParamDisplays[key].html(nf(valueToDisplay, 0, digits));
        }
    }
}

function updateDynamicUIText() {
  select('#stepValueDisplay')?.html(`Step: ${simulationSteps}`);
  select('#timeValueDisplay')?.html(`Time: ${simulationTime.toFixed(2)}`);
  select('#cellCountValueDisplay')?.html(`Cell Count: ${liveCells.size}`);
  select('#maxRadiusValueDisplay')?.html(`Max Radius: ${maxColonyRadius.toFixed(2)}`);

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
  // Click to plant is removed in favor of initial cluster.
  // Could add inspect cell on click later if needed.
  if (mouseX > canvasWidth || mouseY > canvasHeight || mouseX < 0 || mouseY < 0) {
      // Click was outside the main canvas, likely on UI
      return;
  }
   // Example: Log nutrient at clicked cell
    let clickXInGridSpace = mouseX - gridOriginX;
    let clickYInGridSpace = mouseY - gridOriginY;
    const gridCoords = pixelToGridCoords(clickXInGridSpace, clickYInGridSpace); // This needs to be pixelToGridCoords
    const coordStr = `${gridCoords.x},${gridCoords.y}`;

    if (liveCells.has(coordStr)) {
        console.log(`Cell at (${gridCoords.x}, ${gridCoords.y}):`, liveCells.get(coordStr));
        console.log(`Nutrient: ${nutrientGrid.get(coordStr)}`);
    } else {
        console.log(`Empty slot at (${gridCoords.x}, ${gridCoords.y}). Nutrient: ${nutrientGrid.get(coordStr)}`);
    }
}

function keyPressed() { if (key === ' ') { toggleSimulation(); return false; } }

// =============================================================================
