let readstatsSortDirection = {};
let samtoolsSortDirection = {};
let variantsSortDirection = {};
let bedcovSortDirection = {};

function downloadCSV(csvContent, filename) {
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', filename);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

function exportReadstatsToCSV() {
    const table = document.getElementById('readstatsTable');
    if (!table) { alert('Read Statistics table not found.'); return; }
    let csv = [];
    const headers = ['Sample', 'Reads', 'Bases', 'Min_Length', 'Max_Length', 'N50', 'GC_Percent', 'Q20_Percent'];
    csv.push(headers.join(','));
    const dataRows = table.querySelectorAll('tbody tr');
    dataRows.forEach(row => {
        if (row.style.display !== 'none') {
            const cols = [
                row.getAttribute('data-sample'), row.getAttribute('data-reads'), row.getAttribute('data-bases'),
                row.getAttribute('data-minlen'), row.getAttribute('data-maxlen'), row.getAttribute('data-n50'),
                row.getAttribute('data-gc'), row.getAttribute('data-q20')
            ];
            const safeCols = cols.map(text => {
                text = text.trim().replace(/,/g, '');
                if (text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
            });
            csv.push(safeCols.join(','));
        }
    });
    downloadCSV(csv.join('\n'), 'readstats_export.csv');
}

function exportCoverageToCSV() {
    const table = document.getElementById('dataTable');
    let csv = [];
    const headers = ['Sample', 'Chr', 'Gene/Region', 'Region_Size', 'Pct_Ge_1x', 'Pct_Ge_10x', 'Pct_Ge_20x', 'Pct_Ge_30x'];
    csv.push(headers.join(','));
    const dataRows = table.querySelectorAll('tbody tr');
    dataRows.forEach(row => {
        if (row.style.display !== 'none') {
            const cols = [
                row.getAttribute('data-sample'), row.getAttribute('data-chr'), row.getAttribute('data-gene'),
                row.getAttribute('data-total'), row.getAttribute('data-cov1'), row.getAttribute('data-cov10'),
                row.getAttribute('data-cov20'), row.getAttribute('data-cov30')
            ];
            const safeCols = cols.map(text => {
                text = text.trim();
                if (text.includes(',') || text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
            });
            csv.push(safeCols.join(','));
        }
    });
    downloadCSV(csv.join('\n'), 'coverage_breadth_export.csv');
}

function exportSamtoolsToCSV() {
    const table = document.getElementById('samtoolsTable');
    if (!table) { alert('Coverage table not found.'); return; }
    let csv = [];
    const headers = ['Sample', 'Primary_Mapped', 'Primary_Mapped_Pct', 'Bases_On_Target', 'Mean_Target_Coverage', 'Bases_On_Non_Target', 'Mean_Non_Target_Coverage'];
    csv.push(headers.join(','));
    const dataRows = table.querySelectorAll('tbody tr');
    dataRows.forEach(row => {
        if (row.style.display !== 'none') {
            const cols = [
                row.getAttribute('data-sample'),
                row.getAttribute('data-pmapped'), row.getAttribute('data-ppct'),
                row.getAttribute('data-tbases'), row.getAttribute('data-tcov'),
                row.getAttribute('data-ntbases'), row.getAttribute('data-ntcov')
            ];
            const safeCols = cols.map(text => {
                if (!text) return "";
                text = text.trim().replace(/,/g, '');
                if (text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
            });
            csv.push(safeCols.join(','));
        }
    });
    downloadCSV(csv.join('\n'), 'coverage_stats_export.csv');
}

function exportVariantsToCSV() {
    const table = document.getElementById('variantsTable');
    if (!table) { alert('Variants table not found.'); return; }
    let csv = [];
    const headers = ['Sample', 'Total_Variants', 'PASS_Variants', 'SNPs', 'Indels'];
    csv.push(headers.join(','));
    const dataRows = table.querySelectorAll('tbody tr');
    dataRows.forEach(row => {
        if (row.style.display !== 'none') {
            const cols = [
                row.getAttribute('data-sample'),
                row.getAttribute('data-total'), row.getAttribute('data-pass'),
                row.getAttribute('data-snp'), row.getAttribute('data-indel')
            ];
            const safeCols = cols.map(text => {
                if (!text) return "";
                text = text.trim().replace(/,/g, '');
                if (text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
            });
            csv.push(safeCols.join(','));
        }
    });
    downloadCSV(csv.join('\n'), 'variants_export.csv');
}

function exportBedcovToCSV() {
    const table = document.getElementById('bedcovTable');
    if (!table) { alert('Bed Coverage table not found.'); return; }
    let csv = [];
    const headers = ['Sample', 'Chr', 'Gene/Region', 'Start', 'End', 'Size', 'Bases_in_Region', 'Mean_Coverage'];
    csv.push(headers.join(','));
    const dataRows = table.querySelectorAll('tbody tr');
    dataRows.forEach(row => {
        if (row.style.display !== 'none') {
            const cols = [
                row.getAttribute('data-sample'),
                row.getAttribute('data-chr'),
                row.getAttribute('data-gene'),
                row.getAttribute('data-start'),
                row.getAttribute('data-end'),
                row.getAttribute('data-length'),
                row.getAttribute('data-bases'),
                row.getAttribute('data-meancov')
            ];
            const safeCols = cols.map(text => {
                if (!text) return "";
                text = text.trim().replace(/,/g, '');
                if (text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
            });
            csv.push(safeCols.join(','));
        }
    });
    downloadCSV(csv.join('\n'), 'bed_coverage_export.csv');
}

$(document).ready(function () {
    $('#sampleFilter').select2({ placeholder: "Filter samples...", allowClear: true, closeOnSelect: true })
        .on('change', function () {
            filterTable();
            filterReadstatsTable();
            filterSamtoolsTable();
            filterVariantsTable();
            filterBedcovTable();
        });
    $('#regionFilter').select2({ placeholder: "Filter regions...", allowClear: true, closeOnSelect: true })
        .on('change', function () {
            filterTable();
            filterBedcovTable();
        });

    filterTable();
    filterReadstatsTable();
    filterSamtoolsTable();
    filterVariantsTable();
    filterBedcovTable();

    sortDirection[0] = 'desc';
    readstatsSortDirection[0] = 'desc';
    samtoolsSortDirection[0] = 'desc';
    variantsSortDirection[0] = 'desc';
    bedcovSortDirection[0] = 'desc';
    sortTable(0);
    sortReadstatsTable(0);
    sortSamtoolsTable(0);
    sortVariantsTable(0);
    sortBedcovTable(0);
});

function sortReadstatsTable(columnIndex) {
    const table = document.getElementById('readstatsTable');
    if (!table) return;
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.querySelectorAll('tr'));
    readstatsSortDirection[columnIndex] = readstatsSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
    const isAsc = readstatsSortDirection[columnIndex] === 'asc';

    const headers = table.querySelectorAll('th.sortable');
    headers.forEach(h => h.classList.remove('asc', 'desc'));
    if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');

    rows.sort((a, b) => {
        const dataAttrMap = { 0: 'sample', 1: 'reads', 2: 'bases', 4: 'minlen', 5: 'maxlen', 6: 'n50', 7: 'gc', 8: 'q20' };
        const dataKey = dataAttrMap[columnIndex];
        if (!dataKey) {
            const aVal = a.cells[columnIndex].textContent.toLowerCase().trim();
            const bVal = b.cells[columnIndex].textContent.toLowerCase().trim();
            return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        }
        let aVal = a.getAttribute('data-' + dataKey);
        let bVal = b.getAttribute('data-' + dataKey);
        if (columnIndex === 0) return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        return isAsc ? parseFloat(aVal) - parseFloat(bVal) : parseFloat(bVal) - parseFloat(aVal);
    });
    rows.forEach(row => tbody.appendChild(row));
}

function filterReadstatsTable() {
    const selectedSamples = $('#sampleFilter').val() || [];
    const table = document.getElementById('readstatsTable');
    if (!table) return;
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        const sample = row.getAttribute('data-sample');
        const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
        row.style.display = showSample ? '' : 'none';
    });
}

function filterSamtoolsTable() {
    const selectedSamples = $('#sampleFilter').val() || [];
    const table = document.getElementById('samtoolsTable');
    if (!table) return;
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        const sample = row.getAttribute('data-sample');
        const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
        row.style.display = showSample ? '' : 'none';
    });
}

function sortSamtoolsTable(columnIndex) {
    const table = document.getElementById('samtoolsTable');
    if (!table) return;
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.querySelectorAll('tr'));
    samtoolsSortDirection[columnIndex] = samtoolsSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
    const isAsc = samtoolsSortDirection[columnIndex] === 'asc';

    const headers = table.querySelectorAll('th.sortable');
    headers.forEach(h => h.classList.remove('asc', 'desc'));
    if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');

    rows.sort((a, b) => {
        const dataAttrMap = {
            0: 'sample',
            1: 'pmapped', 2: 'ppct',
            3: 'tbases', 4: 'tcov',
            5: 'ntbases', 6: 'ntcov'
        };
        const dataKey = dataAttrMap[columnIndex];
        if (columnIndex === 0) {
            const aVal = a.getAttribute('data-sample');
            const bVal = b.getAttribute('data-sample');
            return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        }
        let aVal = parseFloat(a.getAttribute('data-' + dataKey));
        let bVal = parseFloat(b.getAttribute('data-' + dataKey));
        return isAsc ? aVal - bVal : bVal - aVal;
    });
    rows.forEach(row => tbody.appendChild(row));
}

function sortTable(columnIndex) {
    const table = document.getElementById('dataTable');
    if (!table) return;
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.querySelectorAll('tr'));
    sortDirection[columnIndex] = sortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
    const isAsc = sortDirection[columnIndex] === 'asc';

    const headers = table.querySelectorAll('th.sortable');
    headers.forEach(h => h.classList.remove('asc', 'desc'));
    if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');

    rows.sort((a, b) => {
        const dataAttrMap = { 0: 'sample', 1: 'chr', 2: 'gene', 3: 'total', 4: 'cov1', 5: 'cov10', 6: 'cov20', 7: 'cov30' };
        const dataKey = dataAttrMap[columnIndex];
        let aVal = a.getAttribute('data-' + dataKey);
        let bVal = b.getAttribute('data-' + dataKey);

        if ([0, 1, 2].includes(columnIndex)) {
            return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        } else {
            return isAsc ? parseFloat(aVal) - parseFloat(bVal) : parseFloat(bVal) - parseFloat(aVal);
        }
    });
    rows.forEach(row => tbody.appendChild(row));
}

function filterTable() {
    const selectedSamples = $('#sampleFilter').val() || [];
    const selectedRegions = $('#regionFilter').val() || [];
    const table = document.getElementById('dataTable');
    if (!table) return;
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        const sample = row.getAttribute('data-sample');
        const region = row.getAttribute('data-gene');
        const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
        const showRegion = selectedRegions.length === 0 || selectedRegions.includes(region);
        row.style.display = (showSample && showRegion) ? '' : 'none';
    });
}

function filterVariantsTable() {
    const selectedSamples = $('#sampleFilter').val() || [];
    const table = document.getElementById('variantsTable');
    if (!table) return;
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        const sample = row.getAttribute('data-sample');
        const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
        row.style.display = showSample ? '' : 'none';
    });
}

function sortVariantsTable(columnIndex) {
    const table = document.getElementById('variantsTable');
    if (!table) return;
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.querySelectorAll('tr'));
    variantsSortDirection[columnIndex] = variantsSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
    const isAsc = variantsSortDirection[columnIndex] === 'asc';

    const headers = table.querySelectorAll('th.sortable');
    headers.forEach(h => h.classList.remove('asc', 'desc'));
    if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');

    rows.sort((a, b) => {
        const dataAttrMap = {
            0: 'sample',
            1: 'total', 2: 'pass',
            3: 'snp', 4: 'indel'
        };
        const dataKey = dataAttrMap[columnIndex];
        if (columnIndex === 0) {
            const aVal = a.getAttribute('data-sample');
            const bVal = b.getAttribute('data-sample');
            return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        }
        let aVal = parseFloat(a.getAttribute('data-' + dataKey));
        let bVal = parseFloat(b.getAttribute('data-' + dataKey));
        return isAsc ? aVal - bVal : bVal - aVal;
    });
    rows.forEach(row => tbody.appendChild(row));
}

function filterBedcovTable() {
    const selectedSamples = $('#sampleFilter').val() || [];
    const selectedRegions = $('#regionFilter').val() || [];
    const table = document.getElementById('bedcovTable');
    if (!table) return;
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        const sample = row.getAttribute('data-sample');
        const region = row.getAttribute('data-gene');
        const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
        const showRegion = selectedRegions.length === 0 || selectedRegions.includes(region);
        row.style.display = (showSample && showRegion) ? '' : 'none';
    });
}

function sortBedcovTable(columnIndex) {
    const table = document.getElementById('bedcovTable');
    if (!table) return;
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.querySelectorAll('tr'));
    bedcovSortDirection[columnIndex] = bedcovSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
    const isAsc = bedcovSortDirection[columnIndex] === 'asc';

    const headers = table.querySelectorAll('th.sortable');
    headers.forEach(h => h.classList.remove('asc', 'desc'));
    if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');

    rows.sort((a, b) => {
        const dataAttrMap = {
            0: 'sample',
            1: 'chr',
            2: 'gene',
            3: 'start',
            4: 'end',
            5: 'length',
            6: 'bases',
            7: 'meancov'
        };
        const dataKey = dataAttrMap[columnIndex];

        if ([0, 1, 2].includes(columnIndex)) {
            const aVal = a.getAttribute('data-' + dataKey);
            const bVal = b.getAttribute('data-' + dataKey);
            return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        } else {
            let aVal = parseFloat(a.getAttribute('data-' + dataKey));
            let bVal = parseFloat(b.getAttribute('data-' + dataKey));
            return isAsc ? aVal - bVal : bVal - aVal;
        }
    });
    rows.forEach(row => tbody.appendChild(row));
}
