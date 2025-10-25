import React, { useState } from 'react';
import { LineChart, Line, XAxis, YAxis, Tooltip, CartesianGrid, ResponsiveContainer } from 'recharts';
import { fetchHistoricalData } from './historicalData';
import type { HistoricalHistoryResult } from 'yahoo-finance2/modules/historical';

console.log('App component rendering');

export default function App() {
  const [ticker, setTicker] = useState('');
  const [data, setData] = useState<HistoricalHistoryResult[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const handleLoadData = async () => {
    if (!ticker) return;
    
    setLoading(true);
    setError('');
    
    try {
      const stockData = await fetchHistoricalData(ticker, 30);
      if (!stockData) {
        setError('No data available for this ticker');
      } else {
        setData([stockData]); // Wrap in array since fetchHistoricalData returns a single result
      }
    } catch (e) {
      setError('Error loading data');
      console.error('Error:', e);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{ padding: '20px', backgroundColor: '#f0f0f0', minHeight: '100vh' }}>
      <div style={{ backgroundColor: 'white', padding: '20px', borderRadius: '8px', marginBottom: '20px' }}>
        <h2>Stock Data Viewer</h2>
        <div style={{ display: 'flex', gap: '10px', marginBottom: '10px' }}>
          <input
            type="text"
            value={ticker}
            onChange={(e) => setTicker(e.target.value.toUpperCase())}
            placeholder="Enter stock ticker (e.g. AAPL)"
            style={{ padding: '8px', flexGrow: 1 }}
          />
          <button 
            onClick={handleLoadData}
            disabled={loading || !ticker}
            style={{ padding: '8px 16px', backgroundColor: '#3b82f6', color: 'white', border: 'none', borderRadius: '4px' }}
          >
            {loading ? 'Loading...' : 'Load Data'}
          </button>
        </div>
        {error && <div style={{ color: 'red' }}>{error}</div>}
      </div>

      {data.length > 0 && (
        <div style={{ backgroundColor: 'white', padding: '20px', borderRadius: '8px' }}>
          <h3>{ticker} Stock Price</h3>
          <div style={{ height: '400px' }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="date" />
                <YAxis domain={['auto', 'auto']} />
                <Tooltip />
                <Line type="monotone" dataKey="close" stroke="#3b82f6" dot={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}
    </div>
  );
};
