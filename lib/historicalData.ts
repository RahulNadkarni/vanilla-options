import { HistoricalHistoryResult } from 'yahoo-finance2/modules/historical';

const API_URL = 'http://localhost:3001';

// returns the historical data for a given ticker and date range in days
// if no range is provided, defaults to 30 days
// @param {string} ticker - The stock ticker symbol, e.g. 'AAPL'
// @param {number} range - The date range in days
// @returns {Promise<HistoricalHistoryResult>} - [{ date: Date, open: 150.5, high: 155.2, low: 149.8, close: 154.1, ... }]
export async function fetchHistoricalData(ticker: string, range = 30): Promise<HistoricalHistoryResult> {
    const response = await fetch(`${API_URL}/api/historical/${ticker}?range=${range}`);
    if (!response.ok) {
        throw new Error('Failed to fetch data');
    }
    return response.json();
}