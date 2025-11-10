import { Hono } from 'hono';
import { cors } from 'hono/cors';
import { serve } from '@hono/node-server';
import YahooFinance from 'yahoo-finance2';
import { execFile } from 'node:child_process';
import { promisify } from 'node:util';
import { resolve } from 'node:path';

const app = new Hono();
const yahooFinance = new YahooFinance();
const execFileAsync = promisify(execFile);
const solverBinary = process.env.IV_SOLVER_BIN ?? resolve(process.cwd(), 'build', 'iv_solver_cli');

app.use('*', cors());

app.get('/api/historical/:ticker', async (c) => {
  const ticker = c.req.param('ticker')?.toUpperCase();
  const rangeValue = Number(c.req.query('range'));
  const range = Number.isFinite(rangeValue) && rangeValue > 0 ? rangeValue : 30;

  if (!ticker) {
    return c.json({ error: 'Missing ticker symbol' }, 400);
  }

  try {
    const now = new Date();
    const from = new Date(now.getTime() - range * 86400000);
    const result = await yahooFinance.chart(ticker, {
      period1: from,
      period2: now,
    });

    return c.json({
      meta: result.meta,
      quotes: result.quotes,
    });
  } catch (err: any) {
    return c.json({ error: 'Failed to fetch data', details: err.message ?? String(err) }, 500);
  }
});

app.post('/api/implied-vol', async (c) => {
  try {
    const body = await c.req.json();
    const ticker = typeof body?.ticker === 'string' ? body.ticker.toUpperCase() : '';
    const marketPrice = Number(body?.marketPrice);
    let spot = Number(body?.spot);
    let strike = Number(body?.strike);
    const rate = Number(body?.rate);
    const tenor = Number(body?.tenor);
    let dividend = Number(body?.dividend ?? 0);
    const optionType = typeof body?.optionType === 'string' ? body.optionType.toLowerCase() : '';

    if (!Number.isFinite(dividend)) {
      dividend = 0;
    }

    if (ticker && (!Number.isFinite(spot) || !Number.isFinite(strike))) {
      const now = new Date();
      const from = new Date(now.getTime() - 30 * 86400000);
      const chart = await yahooFinance.chart(ticker, {
        period1: from,
        period2: now,
      });
      const quotes = chart.quotes ?? [];
      const latest = quotes[quotes.length - 1];
      const latestClose = Number(latest?.close ?? latest?.adjclose);
      if (!Number.isFinite(spot) && Number.isFinite(latestClose)) {
        spot = latestClose;
      }
      if (!Number.isFinite(strike) && Number.isFinite(latestClose)) {
        strike = latestClose;
      }
    }

    if (!Number.isFinite(marketPrice) || marketPrice <= 0) {
      return c.json({ error: 'Market price must be positive' }, 400);
    }
    if (!Number.isFinite(spot) || spot <= 0) {
      return c.json({ error: 'Spot must be positive or ticker must have recent data' }, 400);
    }
    if (!Number.isFinite(strike) || strike <= 0) {
      return c.json({ error: 'Strike must be positive or ticker must have recent data' }, 400);
    }
    if (!Number.isFinite(tenor) || tenor <= 0) {
      return c.json({ error: 'Tenor must be positive' }, 400);
    }
    if (!Number.isFinite(rate)) {
      return c.json({ error: 'Rate must be numeric' }, 400);
    }
    if (!Number.isFinite(dividend) || dividend < 0) {
      return c.json({ error: 'Dividend must be zero or positive' }, 400);
    }
    if (optionType !== 'call' && optionType !== 'put') {
      return c.json({ error: 'Invalid option type' }, 400);
    }

    const args = [
      String(marketPrice),
      String(spot),
      String(strike),
      String(rate),
      String(tenor),
      String(dividend),
      optionType,
    ];

    const { stdout } = await execFileAsync(solverBinary, args, { cwd: process.cwd() });
    const parsed = JSON.parse(stdout.trim());
    if (typeof parsed?.impliedVol !== 'number' || !Number.isFinite(parsed.impliedVol)) {
      return c.json({ error: 'Solver returned invalid response' }, 500);
    }
    return c.json({ impliedVol: parsed.impliedVol });
  } catch (err: any) {
    return c.json({ error: 'Implied volatility solver failed', details: err.message ?? String(err) }, 500);
  }
});

const port = Number(process.env.PORT) || 3001;
serve({ fetch: app.fetch, port });
console.log(`Server running on http://localhost:${port}`);
