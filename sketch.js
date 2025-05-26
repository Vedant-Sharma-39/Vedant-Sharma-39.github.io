// --- Core Microbial Strategy Simulation (v32.2 - Perf/Mem Opt & Logic Review) ---
// Model based on "2D Radial Microbial Range Expansion in Fluctuating Environments" description
// Focus: Performance, Memory, Adherence to original model logic.
// ---

// --- Simulation Settings ---
let worldSize = 120;
let cellSize = 4;
let dt = 0.1;
const epsilon = 1e-9;
let drawInterval = 1;

// --- MODEL PARAMETERS (Defaults, adjustable via UI) ---
let PARAMS = {
  W_band: 15,
  initial_cluster_radius: 3,
  default_strategy_type: "Responsive",
  lambda_L_factor: 0.7,
  max_possible_growth_G: 1.0,
  initial_inherent_growth_rate_G: 0.8,
  trade_off_base_lag: 0.5,
  trade_off_max_additional_lag: 3.0,
  trade_off_exponent_k: 2.0,
  lag_mod_factor: 0.5,
  constitutive_growth_G: 0.6,
  constitutive_growth_L: 0.5,
  bet_hedging_p_switch_to_L_if_G_parent: 0.1,
  bet_hedging_p_switch_to_G_if_L_parent: 0.1,
  memory_M_gen_memory_window: 5.0,
  memory_T_lag_GL_multiplier: 0.2,
};

// --- Simulation State ---
let nutrientGrid = new Map();
let liveCells = new Map(); // Key: "x,y", Value: cell object
let nextCellID = 0;
let simulationTime = 0;
let simulationSteps = 0;
let maxColonyRadius = 0;
let previousMaxColonyRadiusForVelocity = 0;
let timeAtPreviousMaxRadiusForVelocity = 0;
let avgRadialVelocity = 0;
let lastStepBirths = 1;
let start = false;
let showNutrientOverlay = false;

// --- Graphics & Grid ---
let nutrientBuffer, cellBuffer, finalStateBuffer;
let finalStateRendered = false;
let canvasWidth, canvasHeight, bufferWidth, bufferHeight, gridOriginX, gridOriginY;

const moore_deltas = [
  { dx: -1, dy: -1 }, { dx: 0, dy: -1 }, { dx: 1, dy: -1 },
  { dx: -1, dy: 0 },  { dx: 1, dy: 0 },
  { dx: -1, dy: 1 },  { dx: 0, dy: 1 },  { dx: 1, dy: 1 },
];

// --- Color Definitions ---
let colorGlucoseBand, colorGalactoseBand;
let colorGSpecialist, colorLSpecialist, colorSwitchingGL;

// --- UI Elements ---
let uiElements = { controlGroups: {}, sliders: {}, selectors: {}, buttons: {} };
let uiParamDisplays = {};
let paramMapping = {};


// =============================================================================
// --- Math & Grid Helpers (Unchanged from v32.1) ---
// =============================================================================
function euclideanDistance(x1, y1, x2, y2) { return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2)); }
function pixelToGridCoords(pixelX, pixelY) { let wx = pixelX - gridOriginX; let wy = pixelY - gridOriginY; return { x: Math.floor((wx + cellSize / 2) / cellSize), y: Math.floor((wy + cellSize / 2) / cellSize) }; }
function gridToPixelCoords(gridX, gridY) { return { x: gridX * cellSize - cellSize / 2, y: gridY * cellSize - cellSize / 2 }; }

// =============================================================================
// --- Drawing Helpers (Unchanged from v32.1) ---
// =============================================================================
function calculateGridDimensions() { bufferWidth = (worldSize * 2 + 1) * cellSize; bufferHeight = (worldSize * 2 + 1) * cellSize; canvasWidth = bufferWidth + 40; canvasHeight = bufferHeight + 40; gridOriginX = bufferWidth / 2; gridOriginY = bufferHeight / 2; }
function drawSquareCell(buffer, gridX, gridY, fillColor, strokeColor = null, weight = 0.5) { const p = gridToPixelCoords(gridX, gridY); buffer.fill(fillColor); if (strokeColor) { buffer.stroke(strokeColor); buffer.strokeWeight(weight); } else { buffer.noStroke(); } buffer.rect(p.x, p.y, cellSize, cellSize); }

// =============================================================================
// --- P5.js Core Functions (Unchanged from v32.1, except minor translation in setup/draw) ---
// =============================================================================
function setup() {
  calculateGridDimensions();
  nutrientBuffer = createGraphics(bufferWidth, bufferHeight);
  cellBuffer = createGraphics(bufferWidth, bufferHeight);
  finalStateBuffer = createGraphics(bufferWidth, bufferHeight);
  finalStateBuffer.pixelDensity(1);

  colorGlucoseBand = color(255, 230, 230, 180); colorGalactoseBand = color(230, 230, 255, 180);
  colorGSpecialist = color(0, 200, 0); colorLSpecialist = color(0, 0, 200); colorSwitchingGL = color(255, 165, 0);

  let totalCanvasHeightUI = Math.max(canvasHeight, 750);
  let mainCanvas = createCanvas(canvasWidth, totalCanvasHeightUI);
  let canvasContainer = select('#canvas-container'); if (canvasContainer) mainCanvas.parent('canvas-container');

  createUI(); initSimulation(); background(240);
  push(); translate( (canvasWidth - bufferWidth) / 2, (canvasHeight - bufferHeight) / 2 );
  image(nutrientBuffer, 0, 0); image(cellBuffer, 0, 0); pop();
  updateUIParamDisplays(); updateDynamicUIText(); noLoop();
}

function draw() {
  if (start) { stepSimulation(); }
  if (frameCount % drawInterval === 0 || !start || showNutrientOverlay) {
    background(240); push(); translate( (canvasWidth - bufferWidth) / 2, (canvasHeight - bufferHeight) / 2 );
    if (showNutrientOverlay) { image(nutrientBuffer, 0, 0); }
    else { image(nutrientBuffer, 0, 0); if (!start && finalStateRendered) { image(finalStateBuffer, 0, 0); } else { image(cellBuffer, 0, 0); } }
    pop();
  }
  updateDynamicUIText();
}

// =============================================================================
// --- Nutrient Environment Logic (Unchanged from v32.1) ---
// =============================================================================
function getNutrientAt(x, y) { const d = euclideanDistance(x,y,0,0); return (Math.floor(d / PARAMS.W_band) % 2 === 0) ? "Glucose" : "Galactose"; }
function populateNutrientGridAndDrawBuffer() { nutrientGrid.clear(); nutrientBuffer.push(); nutrientBuffer.translate(gridOriginX, gridOriginY); nutrientBuffer.background(220); for (let gx = -worldSize; gx <= worldSize; gx++) { for (let gy = -worldSize; gy <= worldSize; gy++) { const nt = getNutrientAt(gx, gy); nutrientGrid.set(`${gx},${gy}`, nt); drawSquareCell(nutrientBuffer, gx, gy, (nt === "Glucose" ? colorGlucoseBand:colorGalactoseBand), null); } } nutrientBuffer.pop(); }

// =============================================================================
// --- Cell Agent Logic (Unchanged from v32.1) ---
// =============================================================================
function trade_off_lag_vs_growth(grG) { const{trade_off_base_lag,trade_off_max_additional_lag,trade_off_exponent_k,max_possible_growth_G}=PARAMS; if(grG<=epsilon)return trade_off_base_lag+trade_off_max_additional_lag; const ngrG=Math.min(grG,max_possible_growth_G)/max_possible_growth_G; return Math.max(0,trade_off_base_lag+Math.pow(ngrG,trade_off_exponent_k)*trade_off_max_additional_lag); }
function calculate_modified_lag(inhLagD,avgInhGparents) { const{lag_mod_factor,max_possible_growth_G}=PARAMS; if(avgInhGparents<epsilon)return inhLagD; const navgF=Math.min(avgInhGparents,max_possible_growth_G)/max_possible_growth_G; return Math.max(0,inhLagD*(1+lag_mod_factor*navgF)); }
function createCell(x,y,strat,grG_val,parentC=null) { nextCellID++; const inhLag_val=trade_off_lag_vs_growth(grG_val); return { id:nextCellID,x:x,y:y,phenotype:"G_specialist",strategy_type:strat,inherent_growth_rate_G:grG_val,inherent_T_lag_GL:inhLag_val,current_growth_rate:0.0,remaining_lag_time:0.0,memory_state:parentC?parentC.memory_state+dt:PARAMS.memory_M_gen_memory_window+1 }; }

// =============================================================================
// --- Initialization & Reset (Unchanged from v32.1) ---
// =============================================================================
function initSimulation() {
  start=false;loop();lastStepBirths=1;finalStateRendered=false;showNutrientOverlay=false; if(uiElements.buttons?.nutrientOverlayBtn)uiElements.buttons.nutrientOverlayBtn.html('Show Nutrients'); if(finalStateBuffer)finalStateBuffer.clear();cellBuffer.clear();
  simulationTime=0;simulationSteps=0;maxColonyRadius=0;nextCellID=0; previousMaxColonyRadiusForVelocity=0;timeAtPreviousMaxRadiusForVelocity=0;avgRadialVelocity=0;
  liveCells.clear(); populateNutrientGridAndDrawBuffer();
  const{initial_cluster_radius,default_strategy_type,initial_inherent_growth_rate_G}=PARAMS; let plantedCount=0;
  for(let gx=-initial_cluster_radius;gx<=initial_cluster_radius;gx++){ for(let gy=-initial_cluster_radius;gy<=initial_cluster_radius;gy++){
    if(euclideanDistance(gx,gy,0,0)<=initial_cluster_radius){ if(nutrientGrid.has(`${gx},${gy}`)){ const coordStr=`${gx},${gy}`; if(!liveCells.has(coordStr)){
      let cell=createCell(gx,gy,default_strategy_type,initial_inherent_growth_rate_G); const local_nutrient=nutrientGrid.get(coordStr);
      if(local_nutrient==="Glucose"){ cell.phenotype="G_specialist"; cell.current_growth_rate=cell.inherent_growth_rate_G; if(cell.strategy_type==="Memory"||(cell.strategy_type==="Responsive"&&PARAMS.memory_M_gen_memory_window>0)) cell.memory_state=PARAMS.memory_M_gen_memory_window+1; }
      else{ cell.phenotype="Switching_GL"; cell.remaining_lag_time=cell.inherent_T_lag_GL; cell.current_growth_rate=0.0; if(cell.strategy_type==="Memory"||(cell.strategy_type==="Responsive"&&PARAMS.memory_M_gen_memory_window>0)) cell.memory_state=0; }
      liveCells.set(coordStr,cell); maxColonyRadius=Math.max(maxColonyRadius,euclideanDistance(gx,gy,0,0)); plantedCount++;
    }}}}
  previousMaxColonyRadiusForVelocity=maxColonyRadius;
  console.log(`Initialized ${plantedCount} cells. Strategy: ${default_strategy_type}. Max Radius: ${maxColonyRadius.toFixed(2)}`);
  drawLiveCellsOnBuffer(cellBuffer); noLoop(); updateDynamicUIText();
}

// =============================================================================
// --- Simulation Step Logic (Optimized Frontier Finding) ---
// =============================================================================
function stepSimulation() {
  if (liveCells.size === 0 || (lastStepBirths === 0 && simulationSteps > 10)) {
    if (start) { let r=liveCells.size===0?"Extinction":"Stasis"; console.log(`${r} @ step ${simulationSteps}, t=${simulationTime.toFixed(2)}. Stopping.`); if(liveCells.size>0&&!finalStateRendered){drawLiveCellsOnBuffer(finalStateBuffer);finalStateRendered=true;} start=false;noLoop();updateDynamicUIText(); } return;
  }
  simulationSteps++; simulationTime += dt; let currentBirthsThisStep = 0;

  // --- 1. Cell State Update (for existing cells) ---
  // Iterate directly over map values, no intermediate array needed if not modifying map size during iteration.
  for (const cell of liveCells.values()) {
    const local_nutrient = nutrientGrid.get(`${cell.x},${cell.y}`);
    let switchedToLagThisStep = false;

    if (cell.phenotype === "Switching_GL") {
      cell.remaining_lag_time -= dt;
      if (cell.remaining_lag_time <= epsilon) {
        cell.remaining_lag_time = 0; cell.phenotype = "L_specialist";
        if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) cell.memory_state = 0;
      }
    } else {
      if (cell.strategy_type === "Responsive" || cell.strategy_type === "BetHedging" || cell.strategy_type === "Memory") {
        if (cell.phenotype === "G_specialist" && local_nutrient === "Galactose") {
          cell.phenotype = "Switching_GL"; switchedToLagThisStep = true;
          if ((cell.strategy_type === "Memory" || cell.strategy_type === "Responsive") && PARAMS.memory_M_gen_memory_window > 0) {
            cell.remaining_lag_time = (cell.memory_state < PARAMS.memory_M_gen_memory_window) ? cell.inherent_T_lag_GL * PARAMS.memory_T_lag_GL_multiplier : cell.inherent_T_lag_GL;
          } else { cell.remaining_lag_time = cell.inherent_T_lag_GL; }
        } else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") { cell.phenotype = "G_specialist"; }
      }
      if (switchedToLagThisStep && (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0))) cell.memory_state = 0;
    }

    cell.current_growth_rate = 0.0;
    if (cell.phenotype === "Switching_GL") { /* Growth 0 */ }
    else if (cell.strategy_type === "Constitutive") cell.current_growth_rate = (local_nutrient === "Glucose") ? PARAMS.constitutive_growth_G : PARAMS.constitutive_growth_L;
    else {
      if (cell.phenotype === "G_specialist" && local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G;
      else if (cell.phenotype === "L_specialist" && local_nutrient === "Galactose") cell.current_growth_rate = PARAMS.lambda_L_factor * cell.inherent_growth_rate_G;
      else if (cell.phenotype === "L_specialist" && local_nutrient === "Glucose") cell.current_growth_rate = cell.inherent_growth_rate_G; // G_specialist on G
    }

    if (cell.strategy_type === "Memory" || (cell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
      if (!((cell.phenotype === "L_specialist" && local_nutrient === "Galactose") || (cell.phenotype === "Switching_GL" && local_nutrient === "Galactose"))) cell.memory_state += dt;
    }
  }

  // --- 2. New Cell Generation (at Frontier) ---
  // Optimized frontier finding and parent assignment
  let frontierSlots = new Map(); // Key: "x,y" of empty slot, Value: { parentalNeighbors: [cellObj,...], sumParentGrowthRate: Z, avgParentInherentGGrowth: X }

  for (const cell of liveCells.values()) { // Iterate over all live cells
    for (const delta of moore_deltas) {   // Check their Moore neighbors
      const nx = cell.x + delta.dx;
      const ny = cell.y + delta.dy;
      // Check bounds first (cheaper)
      if (nx < -worldSize || nx > worldSize || ny < -worldSize || ny > worldSize) continue;
      
      const nCoordStr = `${nx},${ny}`;
      if (!liveCells.has(nCoordStr)) { // If neighbor is an empty slot
        if (!frontierSlots.has(nCoordStr)) {
          frontierSlots.set(nCoordStr, { parentalNeighbors: [], sumParentGrowthRate: 0, sumParentInherentGGrowth: 0, parentCount:0 });
        }
        // The current 'cell' is a parental neighbor to this empty slot 'nCoordStr'
        // This logic was slightly off in 32.1, it collected cell as parent to its own neighbors, which is fine.
        // The model says: "Parental Neighbors: Non-empty cells p_1, p_2, ..., p_k adjacent to (x_slot, y_slot)."
        // So, for each empty slot, we need to find ITS occupied neighbors.
        // The current outer loop (liveCells) finds cells. The inner loop (moore_deltas) finds slots *around* these live cells.
        // This means `cell` is a parent to `nCoordStr`.
        // To correctly implement the model, we should iterate slots, then find parents.
        // OR, the current way: current `cell` IS a parent to `nCoordStr`.
        // The v32.1 version's second pass (iterating frontierSlotsData and then checking *its* neighbors) was more aligned with the model quote.
        // Let's refine:
        // Pass 1: Identify all unique empty frontier slots.
        // Pass 2: For each empty slot, find its adjacent occupied cells (parents).
      }
    }
  }
  
  // Pass 1 (Optimized): Identify unique empty frontier slots by checking neighbors of live cells.
  let emptyFrontierSlotCoords = new Set();
  for (const cell of liveCells.values()) {
      for (const delta of moore_deltas) {
          const nx = cell.x + delta.dx;
          const ny = cell.y + delta.dy;
          if (nx < -worldSize || nx > worldSize || ny < -worldSize || ny > worldSize) continue;
          const nCoordStr = `${nx},${ny}`;
          if (!liveCells.has(nCoordStr)) { // Is an empty slot
              emptyFrontierSlotCoords.add(nCoordStr);
          }
      }
  }

  // Pass 2: For each unique empty slot, find its parental neighbors and their properties.
  for (const slotCoordStr of emptyFrontierSlotCoords) {
      let slotData = { parentalNeighbors: [], sumParentGrowthRate: 0, sumParentInherentGGrowth: 0, parentCount:0 };
      const [sx, sy] = slotCoordStr.split(',').map(Number);

      for (const delta of moore_deltas) { // Check neighbors of the *slot*
          const px = sx + delta.dx;
          const py = sy + delta.dy;
          // Bounds check for parent position (px, py)
          if (px < -worldSize || px > worldSize || py < -worldSize || py > worldSize) continue;

          const pCoordStr = `${px},${py}`;
          if (liveCells.has(pCoordStr)) {
              const parentCell = liveCells.get(pCoordStr);
              slotData.parentalNeighbors.push(parentCell);
              slotData.sumParentGrowthRate += parentCell.current_growth_rate;
              slotData.sumParentInherentGGrowth += parentCell.inherent_growth_rate_G;
              slotData.parentCount++;
          }
      }
      if (slotData.parentCount > 0) {
          frontierSlots.set(slotCoordStr, slotData); // Only add if it has parents
      }
  }


  let newCellsToAdd = []; // {coordStr, cell}
  for (const [slotCoordStr, slotData] of frontierSlots) { // Iterating the populated frontierSlots map
    // Stochastic Colonization: (Unchanged from v32.1, already correct)
    let chosenParent = null;
    if (slotData.sumParentGrowthRate > epsilon) {
      let randVal = random() * slotData.sumParentGrowthRate; let cumulativeGrowth = 0;
      for (const p of slotData.parentalNeighbors) { cumulativeGrowth += p.current_growth_rate; if (randVal <= cumulativeGrowth) { chosenParent = p; break; } }
      if (!chosenParent) chosenParent = slotData.parentalNeighbors[slotData.parentalNeighbors.length -1];
    } else { chosenParent = random(slotData.parentalNeighbors); }
    if (!chosenParent) continue;

    const [nx, ny] = slotCoordStr.split(',').map(Number);
    let daughterCell = createCell(nx, ny, chosenParent.strategy_type, chosenParent.inherent_growth_rate_G, chosenParent);
    const local_nutrient_at_slot = nutrientGrid.get(slotCoordStr);
    daughterCell.phenotype = chosenParent.phenotype;

    if (daughterCell.strategy_type === "BetHedging") { /* ... BetHedging logic (unchanged) ... */ 
        if (chosenParent.phenotype === "G_specialist" && random() < PARAMS.bet_hedging_p_switch_to_L_if_G_parent) daughterCell.phenotype = "L_specialist";
        else if (chosenParent.phenotype === "L_specialist" && random() < PARAMS.bet_hedging_p_switch_to_G_if_L_parent) daughterCell.phenotype = "G_specialist";
    }

    let daughterSwitchedToLag = false;
    if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Galactose") { /* ... Lag Mod logic (unchanged using slotData.sumParentInherentGGrowth / slotData.parentCount) ... */
      daughterCell.phenotype = "Switching_GL"; daughterSwitchedToLag = true;
      let avg_inherent_G_parents = (slotData.parentCount > 0) ? slotData.sumParentInherentGGrowth / slotData.parentCount : 0;
      daughterCell.remaining_lag_time = calculate_modified_lag(daughterCell.inherent_T_lag_GL, avg_inherent_G_parents);
      if ((daughterCell.strategy_type === "Memory" || (daughterCell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) && chosenParent.memory_state < PARAMS.memory_M_gen_memory_window) {
        daughterCell.remaining_lag_time = Math.min(daughterCell.remaining_lag_time, daughterCell.inherent_T_lag_GL * PARAMS.memory_T_lag_GL_multiplier);
      }
    } else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Glucose") { daughterCell.phenotype = "G_specialist"; daughterCell.remaining_lag_time = 0; }
    else { daughterCell.remaining_lag_time = 0; }

    /* ... Daughter's initial growth rate (unchanged) ... */
    daughterCell.current_growth_rate = 0.0;
    if (daughterCell.phenotype === "Switching_GL") { /* Growth 0 */ }
    else if (daughterCell.strategy_type === "Constitutive") daughterCell.current_growth_rate = (local_nutrient_at_slot === "Glucose") ? PARAMS.constitutive_growth_G : PARAMS.constitutive_growth_L;
    else { /* G/L specialist growth */
      if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Glucose") daughterCell.current_growth_rate = daughterCell.inherent_growth_rate_G;
      else if (daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Galactose") daughterCell.current_growth_rate = PARAMS.lambda_L_factor * daughterCell.inherent_growth_rate_G;
      else if (daughterCell.phenotype === "G_specialist" && local_nutrient_at_slot === "Glucose") daughterCell.current_growth_rate = daughterCell.inherent_growth_rate_G;
    }

    /* ... Daughter's memory state init (unchanged) ... */
    if (daughterCell.strategy_type === "Memory" || (daughterCell.strategy_type === "Responsive" && PARAMS.memory_M_gen_memory_window > 0)) {
        if ((daughterCell.phenotype === "L_specialist" && local_nutrient_at_slot === "Galactose") || daughterSwitchedToLag) daughterCell.memory_state = 0;
    }

    newCellsToAdd.push({ coordStr: slotCoordStr, cell: daughterCell });
    currentBirthsThisStep++; maxColonyRadius = Math.max(maxColonyRadius, euclideanDistance(nx, ny, 0, 0));
  }

  newCellsToAdd.forEach(item => liveCells.set(item.coordStr, item.cell));
  lastStepBirths = currentBirthsThisStep;

  if (simulationTime > timeAtPreviousMaxRadiusForVelocity + dt * 10 && maxColonyRadius > previousMaxColonyRadiusForVelocity + 0.1) { /* ... vel calc (unchanged) ... */
      const dr = maxColonyRadius - previousMaxColonyRadiusForVelocity; const dt_vel = simulationTime - timeAtPreviousMaxRadiusForVelocity;
      if (dt_vel > epsilon) avgRadialVelocity = dr / dt_vel;
      previousMaxColonyRadiusForVelocity = maxColonyRadius; timeAtPreviousMaxRadiusForVelocity = simulationTime;
  }
  
  if (start) drawLiveCellsOnBuffer(cellBuffer);
}

// =============================================================================
// --- Drawing Functions (Unchanged from v32.1) ---
// =============================================================================
function drawLiveCellsOnBuffer(targetBuffer) { if(!targetBuffer)return; targetBuffer.clear();targetBuffer.push();targetBuffer.translate(gridOriginX,gridOriginY); for(const[coordStr,cell]of liveCells.entries()){ let cCol; if(cell.phenotype==="Switching_GL")cCol=colorSwitchingGL; else if(cell.phenotype==="G_specialist")cCol=colorGSpecialist; else cCol=colorLSpecialist; drawSquareCell(targetBuffer,cell.x,cell.y,cCol,color(50),0.2); } targetBuffer.pop(); }

// =============================================================================
// --- UI Setup & Updates (Unchanged from v32.1) ---
// =============================================================================
function createUI(){const uiPanel=select('#ui-panel');uiPanel.html('');uiElements={controlGroups:{},sliders:{},selectors:{},buttons:{}};uiParamDisplays={};paramMapping={};const addP=(uK,...pP)=>{paramMapping[uK]=pP;};const getP=(uK)=>{const p=paramMapping[uK];if(!p)return undefined;let v=PARAMS;for(const pp of p){if(v&&v.hasOwnProperty(pp))v=v[pp];else return undefined;}return v;};const setP=(uK,val)=>{const p=paramMapping[uK];if(!p)return;let o=PARAMS;for(let i=0;i<p.length-1;i++){if(o&&o.hasOwnProperty(p[i]))o=o[p[i]];else return;}if(o&&o.hasOwnProperty(p[p.length-1]))o[p[p.length-1]]=val;};const cSR=(lbl,uK,minV,maxV,stV,fmtD=2,...pP)=>{addP(uK,...pP);let currV=getP(uK);if(currV===undefined)return;let rD=createDiv(lbl).parent(uiPanel).addClass('ui-row');let sld=createSlider(minV,maxV,currV,stV).parent(rD).input(()=>{let v=sld.value();setP(uK,v);select('#'+uK+'ValueDisplay').html(nf(v,0,fmtD));if(!start){initSimulation();redraw();}});createSpan(nf(currV,0,fmtD)).parent(rD).id(uK+'ValueDisplay');uiElements.sliders[uK]=sld;uiParamDisplays[uK]=select('#'+uK+'ValueDisplay');};const cSelR=(lbl,uK,opts,...pP)=>{addP(uK,...pP);let currV=getP(uK);if(currV===undefined)return;let rD=createDiv(lbl+": ").parent(uiPanel).addClass('ui-row');let sel=createSelect().parent(rD);opts.forEach(o=>sel.option(typeof o==='string'?o:o.label,typeof o==='string'?o:o.value));sel.selected(currV);sel.changed(()=>{setP(uK,sel.value());if(!start){initSimulation();redraw();}});uiElements.selectors[uK]=sel;};const cSH=(txt)=>createP(`<b>${txt}</b>`).parent(uiPanel);cSH("Simulation Controls:");let ctrlD=createDiv().parent(uiPanel).id('controls-container');uiElements.buttons.startBtn=createButton("Start / Pause (Space)").parent(ctrlD).mousePressed(toggleSimulation);uiElements.buttons.resetBtn=createButton("Reset Sim").parent(ctrlD).mousePressed(()=>{initSimulation();redraw();});uiElements.buttons.nutrientOverlayBtn=createButton("Show Nutrients").parent(ctrlD).mousePressed(toggleNutrientOverlay);cSH("Environment:");cSR("Nutrient Band Width:",'W_band_slider',5,50,1,0,'W_band');cSH("Cell & Strategy Settings:");cSelR("Initial Strategy Type:",'default_strategy_type_sel',["Responsive","Constitutive","BetHedging","Memory"],'default_strategy_type');cSR("Initial Cluster Radius:",'initial_cluster_radius_slider',1,10,1,0,'initial_cluster_radius');cSR("Initial Growth Rate G (lambda_G):",'initial_inherent_growth_rate_G_slider',0.1,1.0,0.01,2,'initial_inherent_growth_rate_G');cSR("Growth Rate L Factor (vs G):",'lambda_L_factor_slider',0.1,1.5,0.01,2,'lambda_L_factor');cSH("G->L Lag Trade-off & Modification:");cSR("Base Lag (T_base):",'trade_off_base_lag_slider',0,2.0,0.05,2,'trade_off_base_lag');cSR("Max Additional Lag (T_add_max):",'trade_off_max_additional_lag_slider',0,5.0,0.1,1,'trade_off_max_additional_lag');cSR("Trade-off Exponent (k):",'trade_off_exponent_k_slider',0.5,4.0,0.1,1,'trade_off_exponent_k');cSR("Parental Fitness Lag Mod Factor:",'lag_mod_factor_slider',-2.0,2.0,0.05,2,'lag_mod_factor');cSH("Strategy: Constitutive");cSR("Constitutive Growth G:",'constitutive_growth_G_slider',0.0,1.0,0.01,2,'constitutive_growth_G');cSR("Constitutive Growth L:",'constitutive_growth_L_slider',0.0,1.0,0.01,2,'constitutive_growth_L');cSH("Strategy: BetHedging");cSR("P(G parent -> L daughter):",'bet_hedging_p_L_slider',0,1,0.01,2,'bet_hedging_p_switch_to_L_if_G_parent');cSR("P(L parent -> G daughter):",'bet_hedging_p_G_slider',0,1,0.01,2,'bet_hedging_p_switch_to_G_if_L_parent');cSH("Strategy: Memory");cSR("Memory Window (time units):",'memory_M_gen_memory_window_slider',0.1,20.0,0.1,1,'memory_M_gen_memory_window');cSR("Memory Lag Multiplier:",'memory_T_lag_GL_multiplier_slider',0.01,1.0,0.01,2,'memory_T_lag_GL_multiplier');cSH("Status:");let statD=createDiv().parent(uiPanel).id('status-section');createSpan('Step: 0').parent(statD).id('stepValueDisplay');createSpan('Time: 0.0').parent(statD).id('timeValueDisplay');createSpan('Cell Count: 0').parent(statD).id('cellCountValueDisplay');createSpan('Max Radius: 0.0').parent(statD).id('maxRadiusValueDisplay');createSpan('Avg. Velocity: 0.00 u/t').parent(statD).id('velocityValueDisplay');createSpan('Status: Initialized').parent(statD).id('statusValueDisplay');updateUIParamDisplays();}
function toggleSimulation(){start=!start;if(start)loop();else noLoop();updateDynamicUIText();if(!start)redraw();}
function toggleNutrientOverlay(){showNutrientOverlay=!showNutrientOverlay;uiElements.buttons.nutrientOverlayBtn.html(showNutrientOverlay?'Show Cells':'Show Nutrients');if(!start)redraw();}
function updateUIParamDisplays(){for(const uK in uiParamDisplays){let valD=paramMapping[uK]?PARAMS[paramMapping[uK][0]]:undefined;if(paramMapping[uK]&¶mMapping[uK].length>1)valD=PARAMS[paramMapping[uK][0]][paramMapping[uK][1]];if(valD===undefined)valD=paramMapping[uK]?getP(paramMapping[uK]):undefined;if(valD===undefined)continue;let digs=2;const pP=paramMapping[uK];const actK=pP[pP.length-1];if(actK.includes('Rate')||actK.includes('Factor')||actK.includes('Multiplier')||actK.includes('Prob')||actK.includes('lag_mod'))digs=2;else if(actK.includes('lag')||actK.includes('exponent')||actK.includes('Window')||actK.includes('max_additional'))digs=1;else if(actK.includes('Radius')||actK.includes('Width')||actK.includes('band'))digs=0;if(uiParamDisplays[uK]&&uiParamDisplays[uK].html()!==nf(valD,0,digs))uiParamDisplays[uK].html(nf(valD,0,digs));}} // simplified getP for flat PARAMS
function getP(path){ if(!path || path.length === 0) return undefined; let val = PARAMS; for(const p of path){ if(val && typeof val === 'object' && val.hasOwnProperty(p)) val = val[p]; else if (PARAMS.hasOwnProperty(p) && path.length === 1) return PARAMS[p]; else return undefined; } return val; } // Refined getP for UI
function updateDynamicUIText(){select('#stepValueDisplay')?.html(`Step: ${simulationSteps}`);select('#timeValueDisplay')?.html(`Time: ${simulationTime.toFixed(2)}`);select('#cellCountValueDisplay')?.html(`Cell Count: ${liveCells.size}`);select('#maxRadiusValueDisplay')?.html(`Max Radius: ${maxColonyRadius.toFixed(2)}`);select('#velocityValueDisplay')?.html(`Avg. Velocity: ${avgRadialVelocity.toFixed(3)} u/t`);let sS=select('#statusValueDisplay');if(sS){let sM="Init",sC="#6c757d";if(start){sM="RUNNING";sC="#28a745";}else if(simulationSteps>0){if(liveCells.size===0){sM="EXTINCT";sC="#dc3545";}else if(lastStepBirths===0&&finalStateRendered){sM="STASIS-Stop";sC="#dc3545";}else if(lastStepBirths===0){sM="PAUSED(Stasis?)";sC="#ffc107";}else{sM="PAUSED";sC="#ffc107";}}else{sM="Ready/PAUSED";sC="#6c757d";}sS.html(`Status: ${sM}`).style('color',sC);}}

// =============================================================================
// --- User Interaction (Unchanged from v32.1) ---
// =============================================================================
function mousePressed(){if(mouseX>(canvasWidth-bufferWidth)/2+bufferWidth||mouseX<(canvasWidth-bufferWidth)/2||mouseY>(canvasHeight-bufferHeight)/2+bufferHeight||mouseY<(canvasHeight-bufferHeight)/2)return;let bX=mouseX-(canvasWidth-bufferWidth)/2;let bY=mouseY-(canvasHeight-bufferHeight)/2;let cXgS=bX-gridOriginX;let cYgS=bY-gridOriginY;const gC=pixelToGridCoords(cXgS+gridOriginX,cYgS+gridOriginY);const cdS=`${gC.x},${gC.y}`;if(liveCells.has(cdS))console.log(`Cell at (${gC.x},${gC.y}):`,liveCells.get(cdS));else console.log(`Empty slot at (${gC.x},${gC.y}).`);console.log(`Nutrient: ${nutrientGrid.get(cdS)}`);}
function keyPressed(){if(key===' ')toggleSimulation();return false;}
// =============================================================================
