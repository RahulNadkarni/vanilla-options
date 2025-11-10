import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  root: 'lib',
  base: '/',
  server: {
    port: 3000,
    open: true
  },
  resolve: {
    extensions: ['.ts', '.tsx', '.js', '.jsx']
  }
})