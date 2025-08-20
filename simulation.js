// --- Cell Types ---
const Empty = 0, Wildtype = 1, Mutant = 2;

// --- Hex Grid Utilities ---
class Hex {
    constructor(q, r, s) { if (Math.round(q + r + s) !== 0) throw new Error("q+r+s=0 not met"); this.q = q; this.r = r; this.s = s; }
    id() { return `${this.q},${this.r}`; }
    neighbors() { return [new Hex(this.q + 1, this.r, this.s - 1), new Hex(this.q - 1, this.r, this.s + 1), new Hex(this.q, this.r + 1, this.s - 1), new Hex(this.q, this.r - 1, this.s + 1), new Hex(this.q + 1, this.r - 1, this.s), new Hex(this.q - 1, this.r + 1, this.s),]; }
}

// --- Plotter (FIXED with Vertical Centering) ---
class HexPlotter {
    constructor(canvas, hexSize, simWidth) {
        this.canvas = canvas; this.ctx = canvas.getContext('2d');
        this.size = hexSize;
        this.colormap = { [Wildtype]: "#003f5c", [Mutant]: "#ff7c43" };
        this.cameraX = 0;

        // FIX: Calculate the vertical offset to center the simulation.
        // The simulation's r-coordinates at q=0 range from approx. 0 to simWidth-1.
        // The midpoint is (simWidth-1)/2. We find its Cartesian y-coordinate.
        const avgR = (simWidth - 1) / 2;
        const centerHexAtQ0 = new Hex(0, avgR, -avgR);
        this.yOffset = this._axialToCartesian(centerHexAtQ0).y;
    }

    _axialToCartesian(h) { return { x: this.size * (3 / 2 * h.q), y: this.size * (Math.sqrt(3) / 2 * h.q + Math.sqrt(3) * h.r) }; }
    _getHexCorners(centerX, centerY) {
        const corners = [];
        for (let i = 0; i < 6; i++) {
            const angleRad = Math.PI / 3 * i;
            corners.push({ x: centerX + this.size * Math.cos(angleRad), y: centerY + this.size * Math.sin(angleRad) });
        }
        return corners;
    }
    _lerp(start, end, amount) { return start + (end - start) * amount; }

    plot(population, targetFrontQ) {
        this.ctx.fillStyle = '#FFFFFF'; this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
        this.ctx.save();
        const { x: targetCameraX } = this._axialToCartesian(new Hex(targetFrontQ, 0, -targetFrontQ));
        this.cameraX = this._lerp(this.cameraX, targetCameraX, 0.08);

        // FIX: Apply the calculated yOffset to the translate function.
        this.ctx.translate(this.canvas.width / 2 - this.cameraX, this.canvas.height / 2 - this.yOffset);

        for (const [id, cellType] of population.entries()) {
            if (cellType === Empty) continue;
            const [q, r] = id.split(',').map(Number);
            const hex = new Hex(q, r, -q - r);
            const { x, y } = this._axialToCartesian(hex);
            const corners = this._getHexCorners(x, y);
            this.ctx.beginPath(); this.ctx.moveTo(corners[0].x, corners[0].y);
            for (let i = 1; i < 6; i++) this.ctx.lineTo(corners[i].x, corners[i].y);
            this.ctx.closePath();
            this.ctx.fillStyle = this.colormap[cellType] || 'gray'; this.ctx.fill();
            this.ctx.strokeStyle = 'white'; this.ctx.lineWidth = 0.5; this.ctx.stroke();
        }
        this.ctx.restore();
    }
}

// --- Gillespie Simulation (Logic is correct, no changes here) ---
class GillespieSimulation {
    constructor(params) {
        this.params = params; this.width = params.width; this.length = params.length;
        this.time = 0; this.stepCount = 0;
        this.mutantCellCount = 0; this.totalCellCount = 0; this.maxQ = 0;
        this.population = new Map(); this.events = []; this.cellToEvents = new Map();
        this._initializePopulation();
        this._findInitialFront();
    }
    updateParams(newParams) { this.params = { ...this.params, ...newParams }; this._rebuildAllEvents(); }
    _initializePopulation() {
        const patchSize = Math.floor(this.width * 0.25);
        const patchStartOffset = Math.floor((this.width - patchSize) / 2);
        for (let i = 0; i < this.width; i++) {
            const rOffset = i;
            const cellType = (i >= patchStartOffset && i < patchStartOffset + patchSize) ? Mutant : Wildtype;
            const h = this._axialToCube(0, rOffset);
            this.population.set(h.id(), cellType);
            if (cellType === Mutant) this.mutantCellCount++;
            this.totalCellCount++;
        }
        this.maxQ = 0;
    }
    _axialToCube(q, rOffset) { const r = rOffset - Math.floor((q + (q & 1)) / 2); return new Hex(q, r, -q - r); }
    _getNeighborsPeriodic(h) {
        return h.neighbors().map(n => {
            const rOffset = n.r + Math.floor((n.q + (n.q & 1)) / 2);
            const wrappedR_Offset = (rOffset % this.width + this.width) % this.width;
            return this._axialToCube(n.q, wrappedR_Offset);
        });
    }
    _getParamsForQ(q) { const { envType, b_m } = this.params; if (envType === 'homogeneous') { return { b_wt: 1.0, b_m }; } else { const patchWidth = 40; const patchIndex = Math.floor(q / patchWidth); if (patchIndex % 2 === 0) return { b_wt: 1.05, b_m: 0.95 }; else return { b_wt: 0.95, b_m: 1.05 }; } }
    _calculateAsymmetricRates() { const { k_total, phi } = this.params; return { k_wt_m: (k_total / 2) * (1 - phi), k_m_wt: (k_total / 2) * (1 + phi) }; }
    _rebuildAllEvents() { this.events = []; this.cellToEvents.clear(); for (const id of this.population.keys()) { this._updateSingleCellEvents(id); } }
    _findInitialFront() { this._rebuildAllEvents(); }
    _updateSingleCellEvents(hexId) {
        const eventIndices = this.cellToEvents.get(hexId) || [];
        for (const index of eventIndices) { this.events[index].rate = 0; }
        this.cellToEvents.delete(hexId);
        const cellType = this.population.get(hexId);
        if (!cellType) return;
        const [q, r] = hexId.split(',').map(Number);
        const h = new Hex(q, r, -q - r);
        const emptyNeighbors = this._getNeighborsPeriodic(h).filter(n => !this.population.has(n.id()) && n.q >= h.q && n.q < this.length);
        if (emptyNeighbors.length > 0) {
            const { k_wt_m, k_m_wt } = this._calculateAsymmetricRates();
            if (cellType === Wildtype && k_wt_m > 0) this._addEvent(k_wt_m, { type: 'switch', parent: h, targetType: Mutant });
            if (cellType === Mutant && k_m_wt > 0) this._addEvent(k_m_wt, { type: 'switch', parent: h, targetType: Wildtype });
            for (const neighbor of emptyNeighbors) {
                const { b_wt, b_m } = this._getParamsForQ(neighbor.q);
                const growthRate = (cellType === Wildtype) ? b_wt : b_m;
                if (growthRate > 0) this._addEvent(growthRate, { type: 'grow', parent: h, target: neighbor });
            }
        }
    }
    _addEvent(rate, data) {
        if (rate <= 0) return; const parentId = data.parent.id();
        if (!this.cellToEvents.has(parentId)) this.cellToEvents.set(parentId, []);
        const eventIndex = this.events.length;
        this.events.push({ rate, ...data }); this.cellToEvents.get(parentId).push(eventIndex);
    }
    _executeEvent(event) {
        if (event.type === 'grow') {
            const parentType = this.population.get(event.parent.id());
            this.population.set(event.target.id(), parentType);
            if (parentType === Mutant) this.mutantCellCount++; this.totalCellCount++;
            if (event.target.q > this.maxQ) this.maxQ = event.target.q;
            this._updateCellAndNeighbors(event.parent); this._updateCellAndNeighbors(event.target);
        } else if (event.type === 'switch') {
            this.population.set(event.parent.id(), event.targetType);
            if (event.targetType === Mutant) this.mutantCellCount++; else this.mutantCellCount--;
            this._updateCellAndNeighbors(event.parent);
        }
    }
    _updateCellAndNeighbors(h) { this._updateSingleCellEvents(h.id()); for (const neighbor of this._getNeighborsPeriodic(h)) { if (this.population.has(neighbor.id())) { this._updateSingleCellEvents(neighbor.id()); } } }
    step() {
        if (this.stepCount % 500 === 0) { this.events = this.events.filter(e => e.rate > 0); this._rebuildCellToEventMap(); }
        const totalRate = this.events.reduce((sum, event) => sum + event.rate, 0);
        if (totalRate <= 1e-9) return false;
        this.time += -Math.log(Math.random()) / totalRate;
        let randVal = Math.random() * totalRate; let chosenEvent = null;
        for (const event of this.events) {
            if (event.rate > 0) { randVal -= event.rate; if (randVal <= 0) { chosenEvent = event; break; } }
        }
        if (chosenEvent) this._executeEvent(chosenEvent);
        this.stepCount++; return this.maxQ < this.length - 2;
    }
    _rebuildCellToEventMap() { this.cellToEvents.clear(); this.events.forEach((event, index) => { const parentId = event.parent.id(); if (!this.cellToEvents.has(parentId)) this.cellToEvents.set(parentId, []); this.cellToEvents.get(parentId).push(index); }); }
}


// --- Main Application Logic ---
document.addEventListener('DOMContentLoaded', () => {
    const canvas = document.getElementById('simulation-canvas');
    canvas.width = 800; canvas.height = 600;
    const controls = { selection: document.getElementById('selection'), k_total_slider: document.getElementById('k_total_slider'), phi: document.getElementById('phi'), envHomogeneous: document.getElementById('env-homogeneous'), envHeterogeneous: document.getElementById('env-heterogeneous'), };
    const valueDisplays = { selection: document.getElementById('selection-value'), k_total: document.getElementById('k_total-value'), phi: document.getElementById('phi-value'), };
    const statDisplays = { time: document.getElementById('time-display'), steps: document.getElementById('steps-display'), cells: document.getElementById('cells-display'), mutantFrac: document.getElementById('mutant-frac-display'), };
    const startPauseBtn = document.getElementById('start-pause');
    const resetBtn = document.getElementById('reset');
    let sim, plotter, isRunning = false, animationFrameId;

    function getKTotal() {
        const slider = controls.k_total_slider; const minLog = -2; const maxLog = 2;
        const position = (slider.value - slider.min) / (slider.max - slider.min);
        const logValue = minLog + position * (maxLog - minLog);
        return Math.pow(10, logValue);
    }
    function getParams() { return { width: 48, length: 2048, b_m: parseFloat(controls.selection.value), k_total: getKTotal(), phi: parseFloat(controls.phi.value), envType: controls.envHomogeneous.checked ? 'homogeneous' : 'heterogeneous', }; }
    function updateValueDisplays() {
        valueDisplays.selection.textContent = parseFloat(controls.selection.value).toFixed(2);
        valueDisplays.phi.textContent = parseFloat(controls.phi.value).toFixed(2);
        const kVal = getKTotal(); valueDisplays.k_total.textContent = kVal < 0.1 ? kVal.toPrecision(2) : kVal.toFixed(2);
    }
    function updateStats() { statDisplays.time.textContent = sim.time.toFixed(2); statDisplays.steps.textContent = sim.stepCount; statDisplays.cells.textContent = sim.totalCellCount; const fraction = sim.totalCellCount > 0 ? sim.mutantCellCount / sim.totalCellCount : 0; statDisplays.mutantFrac.textContent = fraction.toFixed(2); }

    function resetSimulation() {
        isRunning = false; startPauseBtn.textContent = 'Start';
        startPauseBtn.disabled = false; startPauseBtn.classList.remove('running');
        if (animationFrameId) cancelAnimationFrame(animationFrameId);

        sim = new GillespieSimulation(getParams());
        // FIX: Pass the simulation width to the plotter so it can calculate the y-offset.
        plotter = new HexPlotter(canvas, 6, sim.width);

        plotter.plot(sim.population, sim.maxQ);
        updateStats();
    }

    function gameLoop() {
        if (!isRunning) return;
        let continueSim = true; for (let i = 0; i < 10; i++) { if (!sim.step()) { continueSim = false; break; } }
        plotter.plot(sim.population, sim.maxQ); updateStats();
        if (continueSim) { animationFrameId = requestAnimationFrame(gameLoop); }
        else { isRunning = false; startPauseBtn.textContent = 'Finished'; startPauseBtn.disabled = true; }
    }

    startPauseBtn.addEventListener('click', () => { isRunning = !isRunning; if (isRunning) { startPauseBtn.textContent = 'Pause'; startPauseBtn.classList.add('running'); gameLoop(); } else { startPauseBtn.textContent = 'Start'; startPauseBtn.classList.remove('running'); } });
    resetBtn.addEventListener('click', resetSimulation);

    document.querySelectorAll('.controls input').forEach(input => {
        input.addEventListener('input', () => { updateValueDisplays(); if (sim) sim.updateParams(getParams()); });
    });

    updateValueDisplays();
    resetSimulation();
});