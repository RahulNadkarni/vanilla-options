// src/server.ts
import { Hono } from 'hono';
import { cors } from 'hono/cors';
import { serve } from '@hono/node-server';
import YahooFinance from 'yahoo-finance2';

const app = new Hono();
const yahooFinance = new YahooFinance();
app.use('*', cors());

app.get('/api/historical/:ticker', async (c) => {
  const ticker = c.req.param('ticker')?.toUpperCase();
  const range = Number(c.req.query('range')) || 30;

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
    console.error('Fetch error:', err.message);
    return c.json({ error: 'Failed to fetch data', details: err.message }, 500);
  }
});

const port = Number(process.env.PORT) || 3001;
serve({ fetch: app.fetch, port });
console.log(`✅ Server running on http://localhost:${port}`);

// Run with:
// npx hono run src/server.ts
