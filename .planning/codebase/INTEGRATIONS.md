# External Integrations

**Analysis Date:** 2026-02-24

## APIs & External Services

**Data Retrieval:**
- Hosted tree inventory data (treeData.RDS)
  - Service: HTTPS file hosting at doc.ielab.usherbrooke.ca
  - SDK/Client: R base `download.file()` function
  - Used in: `simulations/lambda_plotv0/run_ipm.R` and related simulation scripts
  - URL: https://doc.ielab.usherbrooke.ca/s/qQjXSatrVpGvKzg/download
  - Authentication: Public access (no credentials required)

## Data Storage

**Databases:**
- Not detected - Project does not use external database systems

**File Storage:**
- Local filesystem only
  - Tree data: `treeData.RDS` (R binary serialization format)
  - Parameters: RDS files in `plot_pars/` and `data/output_sim_processed/`
  - Simulation outputs: RDS files in `simulations/*/output/` directories
  - Configuration: `_data.path` file (points to external data directory on local machine)

**Caching:**
- Downloaded treeData.RDS cached locally after first fetch
- SLURM job outputs cached in respective simulation directories

## Authentication & Identity

**Auth Provider:**
- None - Project uses no authentication provider
- Data retrieval is public (HTTPS, no credentials needed)
- No user authentication or identity management

## Monitoring & Observability

**Error Tracking:**
- Not detected

**Logs:**
- SLURM cluster logs (automatically managed by job scheduler)
- R console output captured by SLURM stderr/stdout redirection
- No explicit logging framework implemented

## CI/CD & Deployment

**Hosting:**
- SLURM HPC cluster (job submission via sbatch)
- Local development machine

**CI Pipeline:**
- Not detected - No automated testing, build, or deployment pipeline

**Job Submission:**
- SLURM batch queue system
- Array jobs for parallel simulation execution
- Account: `def-dgravel` (hardcoded in `simulations/lambda_plotv0/run_ipm.R`)
- Job directives: 1 day + 10 hour time limit per job

## Environment Configuration

**Required env vars:**
- BATCH - Batch identifier for job array splitting (custom)
- SLURM_ARRAY_TASK_ID - SLURM automatically assigned array task index

**Configuration files:**
- `_data.path` - Plain text file containing path to external data directory
- No .env files or secrets configuration detected

**Secrets location:**
- No secret management - Project contains no API keys, tokens, or credentials
- External data accessible via public HTTPS link

## Webhooks & Callbacks

**Incoming:**
- None detected

**Outgoing:**
- None detected

## Data Flow

**Initialization:**
1. Download treeData.RDS from remote HTTPS endpoint (one-time download in `run_ipm.R`)
2. Load local RDS files: plot parameters, species-specific posteriors
3. Read simulation control parameters from CSV file

**Execution:**
1. SLURM environment variables set job identifiers (batch_id, array_id)
2. IPM simulation runs on HPC node
3. Output RDS files written to `simulations/*/output/` directory

**Output:**
- Results stored as RDS (R binary format) - no external export detected

---

*Integration audit: 2026-02-24*
