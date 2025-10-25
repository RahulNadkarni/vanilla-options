// src/server.ts
import express, { Request, Response } from 'express';
import cors from 'cors';
import YahooFinance from 'yahoo-finance2';

const app = express();
app.use(cors());
app.use(express.json());

app.get('/api/historical/:ticker', async (req: Request, res: Response) => {
  const ticker = req.params.ticker?.toUpperCase();
  const range = Number(req.query.range) || 30;

  if (!ticker) {
    return res.status(400).json({ error: 'Missing ticker symbol' });
  }

  const yahooFinance = new YahooFinance();

  try {
    const now = new Date();
    const from = new Date(now.getTime() - range * 86400000);

    const result = await yahooFinance.chart(ticker, {
      period1: from,
      period2: now,
    });

    console.log(result)
    res.json(result);
  } catch (err: any) {
    console.error('Fetch error:', err.message);
    res.status(500).json({ error: 'Failed to fetch data', details: err.message });
  }
});

const PORT = Number(process.env.PORT) || 3001;
app.listen(PORT, () => console.log(`✅ Server running on http://localhost:${PORT}`));
