import React, { useState } from 'react';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  Tooltip,
  CartesianGrid,
  ResponsiveContainer,
  Legend,
} from 'recharts';
import { fetchHistoricalData } from './historicalData';
import { solveImpliedVolatility } from './impliedVol';

interface Quote {
  date: string;
  open: number;
  close: number;
  high: number;
  low: number;
  volume: number;
}

const COLORS: Record<string, string> = {
  close: '#3b82f6',
  open: '#10b981',
  high: '#f59e0b',
  low: '#ef4444',
  volume: '#8b5cf6',
};

export default function App() {
  const [ticker, setTicker] = useState('');
  const [range, setRange] = useState(30);
  const [data, setData] = useState<Quote[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [visibleLines, setVisibleLines] = useState({
    close: true,
    open: false,
    high: false,
    low: false,
    volume: false,
  });
  const [ivInputs, setIvInputs] = useState({
    marketPrice: '',
    spot: '',
    strike: '',
    rate: '',
    tenor: '',
    dividend: '',
    optionType: 'call' as 'call' | 'put',
  });
  const [ivLoading, setIvLoading] = useState(false);
  const [ivError, setIvError] = useState('');
  const [ivResult, setIvResult] = useState<number | null>(null);

  const handleLoadData = async () => {
    if (!ticker) return;
    setLoading(true);
    setError('');

    try {
      const result = await fetchHistoricalData(ticker, range);
      const dataObj = result?.quotes || result?.result?.[0]?.quotes;

      if (!dataObj) {
        setError('No data available for this ticker');
        setData([]);
        return;
      }

      const formatted = dataObj.map((q: any) => ({
        ...q,
        date: new Date(q.date).toLocaleDateString(),
      }));

      setData(formatted);
      if (formatted.length > 0) {
        const latest = formatted[formatted.length - 1];
        setIvInputs((prev) => ({
          ...prev,
          spot: latest.close.toString(),
          strike: latest.close.toString(),
        }));
      }
    } catch (e) {
      setError('Error loading data');
      console.error(e);
    } finally {
      setLoading(false);
    }
  };

  const toggleLine = (key: keyof typeof visibleLines) =>
    setVisibleLines((prev) => ({ ...prev, [key]: !prev[key] }));

  const handleIvChange = (key: keyof typeof ivInputs, value: string) => {
    setIvInputs((prev) => ({
      ...prev,
      ...(key === 'optionType'
        ? { optionType: value === 'put' ? 'put' : 'call' }
        : { [key]: value }),
    }));
  };

  const handleSolveIv = async () => {
    setIvError('');
    setIvResult(null);

    const marketPrice = Number(ivInputs.marketPrice);
    const spot = Number(ivInputs.spot);
    const strike = Number(ivInputs.strike);
    const rate = Number(ivInputs.rate);
    const tenor = Number(ivInputs.tenor);
    const dividend = Number(ivInputs.dividend || 0);

    const coreValues = [marketPrice, spot, strike];
    if (!coreValues.every((v) => Number.isFinite(v) && v > 0)) {
      setIvError('Enter valid positive inputs');
      return;
    }
    if (!Number.isFinite(tenor) || tenor <= 0) {
      setIvError('Tenor must be positive');
      return;
    }
    if (!Number.isFinite(rate)) {
      setIvError('Enter a valid rate');
      return;
    }
    if (!Number.isFinite(dividend) || dividend < 0) {
      setIvError('Dividend must be zero or positive');
      return;
    }

    try {
      setIvLoading(true);
      const response = await solveImpliedVolatility({
        marketPrice,
        spot,
        strike,
        rate,
        tenor,
        dividend,
        optionType: ivInputs.optionType,
        ticker: ticker || undefined,
      });
      setIvResult(response.impliedVol);
    } catch (err) {
      setIvError(err instanceof Error ? err.message : 'Failed to solve implied volatility');
    } finally {
      setIvLoading(false);
    }
  };

  const CustomTooltip = ({ active, payload, label }: any) => {
    if (!active || !payload || !payload.length) return null;
    const item = payload[0].payload;
    return (
      <div
        style={{
          background: 'rgba(255,255,255,0.95)',
          border: '1px solid #ddd',
          padding: '8px 12px',
          borderRadius: '8px',
          fontSize: '0.9rem',
        }}
      >
        <div><strong>{label}</strong></div>
        <div>Open: {item.open.toFixed(2)}</div>
        <div>Close: {item.close.toFixed(2)}</div>
        <div>High: {item.high.toFixed(2)}</div>
        <div>Low: {item.low.toFixed(2)}</div>
        <div>Volume: {item.volume.toLocaleString()}</div>
      </div>
    );
  };

  return (
    <div
      style={{
        padding: '30px',
        backgroundColor: '#f3f4f6',
        minHeight: '100vh',
        fontFamily: 'Inter, sans-serif',
      }}
    >
      <div
        style={{
          backgroundColor: 'white',
          padding: '24px',
          borderRadius: '12px',
          boxShadow: '0 2px 6px rgba(0,0,0,0.1)',
          marginBottom: '24px',
        }}
      >
        <h2 style={{ marginBottom: '12px' }}>📈 Stock Data Viewer</h2>

        <div style={{ display: 'flex', gap: '10px', marginBottom: '10px' }}>
          <input
            type="text"
            value={ticker}
            onChange={(e) => setTicker(e.target.value.toUpperCase())}
            placeholder="Enter stock ticker (e.g. AAPL)"
            style={{
              padding: '10px',
              flexGrow: 1,
              borderRadius: '6px',
              border: '1px solid #d1d5db',
            }}
          />
          <input
            type="number"
            value={range}
            onChange={(e) => setRange(Number(e.target.value))}
            min={1}
            max={365}
            style={{
              width: '90px',
              padding: '10px',
              borderRadius: '6px',
              border: '1px solid #d1d5db',
            }}
          />
          <button
            onClick={handleLoadData}
            disabled={loading || !ticker}
            style={{
              padding: '10px 20px',
              backgroundColor: '#3b82f6',
              color: 'white',
              border: 'none',
              borderRadius: '6px',
              fontWeight: 500,
              cursor: 'pointer',
            }}
          >
            {loading ? 'Loading...' : 'Load Data'}
          </button>
        </div>

        <div style={{ display: 'flex', gap: '10px', flexWrap: 'wrap' }}>
          {Object.keys(visibleLines).map((key) => (
            <button
              key={key}
              onClick={() => toggleLine(key as keyof typeof visibleLines)}
              style={{
                padding: '6px 12px',
                borderRadius: '6px',
                border: visibleLines[key as keyof typeof visibleLines]
                  ? `2px solid ${COLORS[key]}`
                  : '1px solid #d1d5db',
                backgroundColor: visibleLines[key as keyof typeof visibleLines]
                  ? COLORS[key] + '22'
                  : 'white',
                cursor: 'pointer',
                fontSize: '0.85rem',
              }}
            >
              {visibleLines[key as keyof typeof visibleLines]
                ? `✅ ${key}`
                : key}
            </button>
          ))}
        </div>

        {error && <div style={{ color: 'red', marginTop: '10px' }}>{error}</div>}
      </div>

      <div
        style={{
          backgroundColor: 'white',
          padding: '24px',
          borderRadius: '12px',
          boxShadow: '0 2px 6px rgba(0,0,0,0.1)',
          marginBottom: '24px',
        }}
      >
        <h2 style={{ marginBottom: '12px' }}>🧮 Implied Volatility Solver</h2>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(160px, 1fr))', gap: '12px' }}>
          <input
            type="number"
            placeholder="Market Price"
            value={ivInputs.marketPrice}
            onChange={(e) => handleIvChange('marketPrice', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          />
          <input
            type="number"
            placeholder="Spot"
            value={ivInputs.spot}
            onChange={(e) => handleIvChange('spot', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          />
          <input
            type="number"
            placeholder="Strike"
            value={ivInputs.strike}
            onChange={(e) => handleIvChange('strike', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          />
          <input
            type="number"
            placeholder="Rate"
            value={ivInputs.rate}
            onChange={(e) => handleIvChange('rate', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          />
          <input
            type="number"
            placeholder="Tenor (years)"
            value={ivInputs.tenor}
            onChange={(e) => handleIvChange('tenor', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          />
          <input
            type="number"
            placeholder="Dividend"
            value={ivInputs.dividend}
            onChange={(e) => handleIvChange('dividend', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          />
          <select
            value={ivInputs.optionType}
            onChange={(e) => handleIvChange('optionType', e.target.value)}
            style={{ padding: '10px', borderRadius: '6px', border: '1px solid #d1d5db' }}
          >
            <option value="call">Call</option>
            <option value="put">Put</option>
          </select>
        </div>
        <button
          onClick={handleSolveIv}
          disabled={ivLoading}
          style={{
            padding: '10px 20px',
            backgroundColor: '#0f766e',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            fontWeight: 500,
            cursor: 'pointer',
            marginTop: '16px',
          }}
        >
          {ivLoading ? 'Solving...' : 'Solve'}
        </button>
        {ivError && <div style={{ color: 'red', marginTop: '10px' }}>{ivError}</div>}
        {ivResult !== null && (
          <div style={{ marginTop: '12px', fontSize: '1.1rem' }}>
            Implied volatility: <strong>{ivResult.toFixed(6)}</strong>
          </div>
        )}
      </div>

      {data.length > 0 && (
        <div
          style={{
            backgroundColor: 'white',
            padding: '24px',
            borderRadius: '12px',
            boxShadow: '0 2px 6px rgba(0,0,0,0.1)',
          }}
        >
          <h3 style={{ marginBottom: '12px' }}>
            {ticker} — Last {range} Days
          </h3>
          <div style={{ height: '450px' }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="date" />
                <YAxis />
                <Tooltip content={<CustomTooltip />} />
                <Legend />
                {Object.entries(visibleLines).map(
                  ([key, visible]) =>
                    visible && (
                      <Line
                        key={key}
                        type="monotone"
                        dataKey={key}
                        stroke={COLORS[key]}
                        dot={false}
                        strokeWidth={2}
                      />
                    )
                )}
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}
    </div>
  );
}
