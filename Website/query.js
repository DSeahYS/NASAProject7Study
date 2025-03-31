class DataQuery {
    constructor() {
        this.form = document.getElementById('query-form');
        this.resultsContent = document.getElementById('results-content');
        this.validationContent = document.getElementById('validation-content');
        this.clearBtn = document.getElementById('clear-results');
        this.tabBtns = document.querySelectorAll('.tab-btn');
        this.tabContents = document.querySelectorAll('.tab-content');
        
        this.setupEventListeners();
    }
    
    setupEventListeners() {
        this.form.addEventListener('submit', (e) => this.handleSubmit(e));
        this.clearBtn.addEventListener('click', () => this.clearResults());
        
        this.tabBtns.forEach(btn => {
            btn.addEventListener('click', () => this.switchTab(btn));
        });
    }
    
    async handleSubmit(e) {
        e.preventDefault();
        
        const formData = new FormData(this.form);
        const params = {
            objectId: formData.get('object-id'),
            radius: formData.get('search-radius'),
            instruments: Array.from(document.querySelectorAll('input[name="instrument"]:checked'))
                .map(el => el.value),
            bands: Array.from(document.querySelectorAll('input[name="band"]:checked'))
                .map(el => el.value),
            validate: document.getElementById('validate-data').checked
        };
        
        if (!params.objectId) {
            alert('Please enter an object ID or coordinates');
            return;
        }
        
        try {
            // In a real app, this would call your API
            // const response = await fetch('/api/query', {
            //     method: 'POST',
            //     headers: { 'Content-Type': 'application/json' },
            //     body: JSON.stringify(params)
            // });
            // const data = await response.json();
            
            // Mock response for demonstration
            const data = this.generateMockResponse(params);
            
            this.displayResults(data.results);
            
            if (params.validate && data.validation) {
                this.displayValidation(data.validation);
                document.querySelector('.tab-btn[data-tab="validation"]').click();
            }
            
        } catch (error) {
            this.showError('Query failed', error);
        }
    }
    
    generateMockResponse(params) {
        // Generate mock data based on query parameters
        const results = [];
        const numResults = Math.floor(Math.random() * 5) + 1;
        
        for (let i = 0; i < numResults; i++) {
            results.push({
                instrument: params.instruments[Math.floor(Math.random() * params.instruments.length)],
                band: params.bands[Math.floor(Math.random() * params.bands.length)],
                date: new Date().toISOString(),
                exposure: (Math.random() * 1000 + 100).toFixed(0) + ' seconds',
                url: `https://archive.example.com/data/${params.objectId}_${i}.fits`
            });
        }
        
        const response = { results };
        
        if (params.validate) {
            response.validation = {
                astrometry: {
                    status: 'Success',
                    total_rms: (Math.random() * 0.5 + 0.1).toFixed(2) + ' arcsec',
                    offset: (Math.random() * 0.3 + 0.05).toFixed(2) + ' arcsec'
                },
                photometry: {
                    status: 'Success',
                    zp_offset: (Math.random() * 0.2 + 0.01).toFixed(2) + ' mag',
                    mag_offset: (Math.random() * 0.15 + 0.01).toFixed(2) + ' mag'
                },
                matching: {
                    status: 'Success',
                    n_matches: Math.floor(Math.random() * 50 + 10)
                }
            };
        }
        
        return response;
    }
    
    displayResults(results) {
        if (results.length === 0) {
            this.resultsContent.innerHTML = '<p>No results found for this query.</p>';
            return;
        }
        
        let html = '<div class="results-grid">';
        
        results.forEach((result, i) => {
            html += `
                <div class="result-item">
                    <h3>Result #${i + 1}</h3>
                    <p><strong>Instrument:</strong> ${result.instrument}</p>
                    <p><strong>Band:</strong> ${result.band}</p>
                    <p><strong>Date:</strong> ${new Date(result.date).toLocaleString()}</p>
                    <p><strong>Exposure:</strong> ${result.exposure}</p>
                    <p><strong>URL:</strong> <a href="${result.url}" target="_blank">${result.url}</a></p>
                </div>
                <hr>
            `;
        });
        
        html += '</div>';
        this.resultsContent.innerHTML = html;
    }
    
    displayValidation(validation) {
        let html = `
            <h3>GAIA DR3 Validation Report</h3>
            <div class="validation-section">
                <h4>Astrometry</h4>
                <p><strong>Status:</strong> ${validation.astrometry.status}</p>
                <p><strong>Total RMS:</strong> ${validation.astrometry.total_rms}</p>
                <p><strong>Positional Offset:</strong> ${validation.astrometry.offset}</p>
            </div>
            <div class="validation-section">
                <h4>Photometry</h4>
                <p><strong>Status:</strong> ${validation.photometry.status}</p>
                <p><strong>Zero Point Offset:</strong> ${validation.photometry.zp_offset}</p>
                <p><strong>Magnitude Offset:</strong> ${validation.photometry.mag_offset}</p>
            </div>
            <div class="validation-section">
                <h4>Source Matching</h4>
                <p><strong>Status:</strong> ${validation.matching.status}</p>
                <p><strong>Matched Sources:</strong> ${validation.matching.n_matches}</p>
            </div>
        `;
        
        this.validationContent.innerHTML = html;
    }
    
    clearResults() {
        this.resultsContent.innerHTML = '<p>No query results yet. Enter parameters and click "Execute Query".</p>';
        this.validationContent.innerHTML = '<p>No validation results available. Enable validation and run a query.</p>';
    }
    
    switchTab(btn) {
        const tabId = btn.getAttribute('data-tab');
        
        this.tabBtns.forEach(b => b.classList.remove('active'));
        this.tabContents.forEach(c => c.classList.remove('active'));
        
        btn.classList.add('active');
        document.getElementById(`${tabId}-tab`).classList.add('active');
    }
    
    showError(message, error) {
        console.error(message, error);
        this.resultsContent.innerHTML = `
            <div class="error-message">
                <h3>${message}</h3>
                <p>${error.message}</p>
            </div>
        `;
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new DataQuery();
});