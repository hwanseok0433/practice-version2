    // Initialize Lucide Icons
    lucide.createIcons();

    // Theme toggler
    function toggleDarkMode() {
      const isDark = document.documentElement.classList.toggle('dark');
      document.getElementById('sun-icon').classList.toggle('hidden', !isDark);
      document.getElementById('moon-icon').classList.toggle('hidden', isDark);
      // Update charts for theme
      updateSimulationChartsTheme();
    }

    // Tab Navigation
    function switchTab(tabId) {
      // Hide all tabs
      document.querySelectorAll('.tab-content').forEach(el => {
        el.classList.add('hidden');
      });

      // Show selected tab
      const target = document.getElementById('tab-' + tabId);
      if (target) {
        target.classList.remove('hidden');
      }

      // Update button visual styles
      document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.classList.remove('bg-white', 'dark:bg-slate-700', 'text-blue-600', 'dark:text-blue-400', 'shadow-sm');
        btn.classList.add('text-slate-600', 'dark:text-slate-400');
      });

      const activeBtn = document.getElementById('nav-' + tabId);
      if (activeBtn) {
        activeBtn.classList.add('bg-white', 'dark:bg-slate-700', 'text-blue-600', 'dark:text-blue-400', 'shadow-sm');
        activeBtn.classList.remove('text-slate-600', 'dark:text-slate-400');
      }

      // Tab specific re-render/resize
      if (tabId === 'simulation') {
        setTimeout(() => {
          try {
            updateSimulation();
            if (chartParallel) chartParallel.resize();
            if (chartCounter) chartCounter.resize();
          } catch (e) {
            console.warn(e);
          }
        }, 50);
      } else if (tabId === 'lab') {
        setTimeout(() => {
          try {
            calcLabResults();
            if (chartLabParallel) chartLabParallel.resize();
            if (chartLabCounter) chartLabCounter.resize();
          } catch (e) {
            console.warn(e);
          }
        }, 50);
      }

      if (window.lucide && lucide.createIcons) {
        lucide.createIcons();
      }
    }

    // ==========================================
    // PHYSICAL CONSTANTS (Slide 8 & 9)
    // ==========================================
    const CONSTANTS = {
      Cp: 4184, // J/(kg*K) = cal/(g*C) * 4.184
      rho: 1.0, // kg/L
      Di: 0.01426, // m (14.26 mm)
      Do: 0.01588, // m (15.88 mm)
      L_default: 1.5, // m (1500 mm)
      Ai: Math.PI * 0.01426 * 1.5, // 0.067198 m2
      Ao: Math.PI * 0.01588 * 1.5, // 0.07483 m2
      Am: (Math.PI * (0.01588 - 0.01426) * 1.5) / Math.log(0.01588 / 0.01426) // 0.0709 m2
    };

    // ==========================================
    // SIMULATION ENGINE (Analytical & NTU method)
    // ==========================================
    let chartParallel = null;
    let chartCounter = null;
    let chartLabParallel = null;
    let chartLabCounter = null;

    function solveHeatExchanger(Th_in, Tc_in, mh_Lmin, mc_Lmin, U, L) {
      // Flow rates: L/min -> kg/s
      const mh = (mh_Lmin * CONSTANTS.rho) / 60.0;
      const mc = (mc_Lmin * CONSTANTS.rho) / 60.0;
      
      const Ch = mh * CONSTANTS.Cp; // W/K
      const Cc = mc * CONSTANTS.Cp; // W/K

      const Cmin = Math.min(Ch, Cc);
      const Cmax = Math.max(Ch, Cc);
      const Cr = Cmin / Cmax;

      const A = Math.PI * CONSTANTS.Do * L; // outer heat transfer area
      const NTU = (U * A) / Cmin;

      // 1. Parallel Flow Effectiveness
      const eps_par = (1 - Math.exp(-NTU * (1 + Cr))) / (1 + Cr);
      const q_max = Cmin * (Th_in - Tc_in);
      const q_par = eps_par * q_max;

      const Th_out_par = Th_in - q_par / Ch;
      const Tc_out_par = Tc_in + q_par / Cc;

      const dt1_par = Th_in - Tc_in;
      const dt2_par = Th_out_par - Tc_out_par;
      const lmtd_par = (dt1_par === dt2_par) ? dt1_par : (dt1_par - dt2_par) / Math.log(dt1_par / dt2_par);

      // 2. Counterflow Effectiveness
      let eps_cnt;
      if (Cr < 0.9999) {
        const expTerm = Math.exp(-NTU * (1 - Cr));
        eps_cnt = (1 - expTerm) / (1 - Cr * expTerm);
      } else {
        eps_cnt = NTU / (1 + NTU);
      }
      const q_cnt = eps_cnt * q_max;

      const Th_out_cnt = Th_in - q_cnt / Ch;
      const Tc_out_cnt = Tc_in + q_cnt / Cc;

      const dt1_cnt = Th_in - Tc_out_cnt;
      const dt2_cnt = Th_out_cnt - Tc_in;
      const lmtd_cnt = (Math.abs(dt1_cnt - dt2_cnt) < 1e-4) ? dt1_cnt : (dt1_cnt - dt2_cnt) / Math.log(dt1_cnt / dt2_cnt);

      // Profile along length: N steps
      const N = 20;
      const xVals = [];
      const par_Th_curve = [];
      const par_Tc_curve = [];
      const cnt_Th_curve = [];
      const cnt_Tc_curve = [];

      for (let i = 0; i <= N; i++) {
        const x = (i / N) * L;
        xVals.push(x.toFixed(2));

        // Local UA up to x
        const Ax = Math.PI * CONSTANTS.Do * x;

        // Parallel local profile
        // Th(x) - Tc(x) = (Th_in - Tc_in) * exp(-U * Ax * (1/Ch + 1/Cc))
        const decayPar = Math.exp(-U * Ax * (1/Ch + 1/Cc));
        const deltaT_par_x = (Th_in - Tc_in) * decayPar;
        const q_par_x = ( (Th_in - Tc_in) - deltaT_par_x ) / (1/Ch + 1/Cc);
        par_Th_curve.push((Th_in - q_par_x / Ch).toFixed(1));
        par_Tc_curve.push((Tc_in + q_par_x / Cc).toFixed(1));

        // Counterflow local profile
        // Tc(x) and Th(x)
        if (Cr < 0.9999) {
          const exponent = U * Ax * (1/Cc - 1/Ch);
          const deltaT_cnt_x = (Th_out_cnt - Tc_in) * Math.exp(exponent);
          const Th_x = Th_out_cnt + ( (deltaT_cnt_x - (Th_out_cnt - Tc_in)) / (1 - Ch/Cc) );
          const Tc_x = Tc_in + (Ch/Cc) * (Th_x - Th_out_cnt);
          cnt_Th_curve.push(Th_x.toFixed(1));
          cnt_Tc_curve.push(Tc_x.toFixed(1));
        } else {
          // Cr = 1 linear profile
          const deltaT_cnt = Th_in - Tc_out_cnt;
          const q_rate = (U * A * deltaT_cnt) / L;
          const Th_x = Th_in - (q_rate * x) / Ch;
          const Tc_x = Tc_out_cnt - (q_rate * x) / Cc;
          cnt_Th_curve.push(Th_x.toFixed(1));
          cnt_Tc_curve.push(Tc_x.toFixed(1));
        }
      }

      return {
        parallel: {
          Th_out: Th_out_par,
          Tc_out: Tc_out_par,
          q: q_par,
          lmtd: lmtd_par,
          eps: eps_par,
          Th_curve: par_Th_curve,
          Tc_curve: par_Tc_curve
        },
        counter: {
          Th_out: Th_out_cnt,
          Tc_out: Tc_out_cnt,
          q: q_cnt,
          lmtd: lmtd_cnt,
          eps: eps_cnt,
          Th_curve: cnt_Th_curve,
          Tc_curve: cnt_Tc_curve
        },
        xVals: xVals
      };
    }

    function updateSimulation() {
      const Th_in = parseFloat(document.getElementById('sim-Th-in').value);
      const Tc_in = parseFloat(document.getElementById('sim-Tc-in').value);
      const mh = parseFloat(document.getElementById('sim-mh').value);
      const mc = parseFloat(document.getElementById('sim-mc').value);
      const U = parseFloat(document.getElementById('sim-U-coeff').value);
      const L_mm = parseFloat(document.getElementById('sim-L-pipe').value);
      const L = L_mm / 1000.0;

      // Update badge labels
      document.getElementById('val-Th-in').innerText = Th_in.toFixed(1) + ' ℃';
      document.getElementById('val-Tc-in').innerText = Tc_in.toFixed(1) + ' ℃';
      document.getElementById('val-mh').innerText = mh.toFixed(1) + ' L/min';
      document.getElementById('val-mc').innerText = mc.toFixed(1) + ' L/min';
      document.getElementById('val-U-coeff').innerText = U.toLocaleString() + ' W/m²·K';
      document.getElementById('val-L-pipe').innerText = L_mm.toLocaleString() + ' mm';

      const res = solveHeatExchanger(Th_in, Tc_in, mh, mc, U, L);

      // Update Parallel numbers
      document.getElementById('par-Th-out').innerText = res.parallel.Th_out.toFixed(1) + ' ℃';
      document.getElementById('par-Tc-out').innerText = res.parallel.Tc_out.toFixed(1) + ' ℃';
      document.getElementById('par-q').innerText = Math.round(res.parallel.q).toLocaleString() + ' W';
      document.getElementById('par-lmtd').innerText = res.parallel.lmtd.toFixed(1) + ' ℃';

      // Update Counterflow numbers
      document.getElementById('cnt-Th-out').innerText = res.counter.Th_out.toFixed(1) + ' ℃';
      document.getElementById('cnt-Tc-out').innerText = res.counter.Tc_out.toFixed(1) + ' ℃';
      document.getElementById('cnt-q').innerText = Math.round(res.counter.q).toLocaleString() + ' W';
      document.getElementById('cnt-lmtd').innerText = res.counter.lmtd.toFixed(1) + ' ℃';

      // Advantage comparison
      const dutyDiff = ((res.counter.q - res.parallel.q) / res.parallel.q) * 100;
      document.getElementById('comp-duty-diff').innerText = `+${dutyDiff.toFixed(1)}% 더 많은 열전달`;

      // Check temperature crossover
      const hasCross = res.counter.Tc_out > res.counter.Th_out;
      document.getElementById('cross-badge').classList.toggle('hidden', !hasCross);

      // Render Charts
      renderParallelChart(res.xVals, res.parallel.Th_curve, res.parallel.Tc_curve);
      renderCounterChart(res.xVals, res.counter.Th_curve, res.counter.Tc_curve);
    }

    function renderParallelChart(labels, hotData, coldData) {
      if (!window.Chart) return;
      const canvas = document.getElementById('chart-parallel');
      if (!canvas) return;
      const ctx = canvas.getContext('2d');
      if (chartParallel) {
        chartParallel.data.labels = labels;
        chartParallel.data.datasets[0].data = hotData;
        chartParallel.data.datasets[1].data = coldData;
        chartParallel.update('none');
        return;
      }
      chartParallel = new Chart(ctx, {
        type: 'line',
        data: {
          labels: labels,
          datasets: [
            {
              label: '온수 (고온유체)',
              data: hotData,
              borderColor: '#f97316', // Orange
              backgroundColor: 'rgba(249, 115, 22, 0.1)',
              borderWidth: 3,
              tension: 0.35,
              fill: false,
              pointRadius: 2
            },
            {
              label: '냉각수 (저온유체)',
              data: coldData,
              borderColor: '#3b82f6', // Blue
              backgroundColor: 'rgba(59, 130, 246, 0.1)',
              borderWidth: 3,
              tension: 0.35,
              fill: false,
              pointRadius: 2
            }
          ]
        },
        options: getChartOptions('관 길이 x (m)', '온도 (℃)')
      });
    }

    function renderCounterChart(labels, hotData, coldData) {
      if (!window.Chart) return;
      const canvas = document.getElementById('chart-counter');
      if (!canvas) return;
      const ctx = canvas.getContext('2d');
      if (chartCounter) {
        chartCounter.data.labels = labels;
        chartCounter.data.datasets[0].data = hotData;
        chartCounter.data.datasets[1].data = coldData;
        chartCounter.update('none');
        return;
      }
      chartCounter = new Chart(ctx, {
        type: 'line',
        data: {
          labels: labels,
          datasets: [
            {
              label: '온수 (고온유체)',
              data: hotData,
              borderColor: '#f97316', // Orange
              backgroundColor: 'rgba(249, 115, 22, 0.1)',
              borderWidth: 3,
              tension: 0.35,
              fill: false,
              pointRadius: 2
            },
            {
              label: '냉각수 (저온유체)',
              data: coldData,
              borderColor: '#3b82f6', // Blue
              backgroundColor: 'rgba(59, 130, 246, 0.1)',
              borderWidth: 3,
              tension: 0.35,
              fill: false,
              pointRadius: 2
            }
          ]
        },
        options: getChartOptions('관 길이 x (m)', '온도 (℃)')
      });
    }

    function getChartOptions(xTitle, yTitle) {
      const isDark = document.documentElement.classList.contains('dark');
      const textColor = isDark ? '#94a3b8' : '#64748b';
      const gridColor = isDark ? '#334155' : '#f1f5f9';

      return {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: { display: false },
          tooltip: {
            mode: 'index',
            intersect: false,
            callbacks: {
              label: function(context) {
                return `${context.dataset.label}: ${context.parsed.y} ℃`;
              }
            }
          }
        },
        scales: {
          x: {
            title: { display: true, text: xTitle, color: textColor, font: { size: 10 } },
            ticks: { color: textColor, font: { size: 10 } },
            grid: { color: gridColor }
          },
          y: {
            title: { display: true, text: yTitle, color: textColor, font: { size: 10 } },
            ticks: { color: textColor, font: { size: 10 } },
            grid: { color: gridColor }
          }
        }
      };
    }

    function updateSimulationChartsTheme() {
      if (chartParallel) {
        chartParallel.options = getChartOptions('관 길이 x (m)', '온도 (℃)');
        chartParallel.update();
      }
      if (chartCounter) {
        chartCounter.options = getChartOptions('관 길이 x (m)', '온도 (℃)');
        chartCounter.update();
      }
      if (chartLabParallel) {
        chartLabParallel.options = getChartOptions('위치 (입구 → 출구)', '온도 (℃)');
        chartLabParallel.update();
      }
      if (chartLabCounter) {
        chartLabCounter.options = getChartOptions('위치 (입구 → 출구)', '온도 (℃)');
        chartLabCounter.update();
      }
    }

    // Presets for Simulation
    function applySimPreset(type) {
      if (type === 'standard') {
        document.getElementById('sim-Th-in').value = 58.0;
        document.getElementById('sim-Tc-in').value = 15.0;
        document.getElementById('sim-mh').value = 2.0;
        document.getElementById('sim-mc').value = 3.0;
        document.getElementById('sim-U-coeff').value = 1800;
        document.getElementById('sim-L-pipe').value = 1500;
      } else if (type === 'highHotFlow') {
        document.getElementById('sim-Th-in').value = 65.0;
        document.getElementById('sim-Tc-in').value = 12.0;
        document.getElementById('sim-mh').value = 6.5;
        document.getElementById('sim-mc').value = 2.0;
        document.getElementById('sim-U-coeff').value = 2200;
        document.getElementById('sim-L-pipe').value = 1500;
      } else if (type === 'tempCross') {
        document.getElementById('sim-Th-in').value = 75.0;
        document.getElementById('sim-Tc-in').value = 20.0;
        document.getElementById('sim-mh').value = 1.2;
        document.getElementById('sim-mc').value = 1.0;
        document.getElementById('sim-U-coeff').value = 2400;
        document.getElementById('sim-L-pipe').value = 2500;
      }
      updateSimulation();
    }

    // ==========================================
    // LAB DATA CALCULATION ENGINE (Slide 9 & 10)
    // ==========================================
    function getAverage(ids) {
      const vals = ids.map(id => parseFloat(document.getElementById(id).value) || 0);
      const sum = vals.reduce((a, b) => a + b, 0);
      return sum / vals.length;
    }

    function calcLabResults() {
      // 1. Counterflow Averages
      const cnt_T1 = getAverage(['cnt-t1-1', 'cnt-t1-2', 'cnt-t1-3']);
      const cnt_T2 = getAverage(['cnt-t2-1', 'cnt-t2-2', 'cnt-t2-3']);
      const cnt_t1 = getAverage(['cnt-tc1-1', 'cnt-tc1-2', 'cnt-tc1-3']); // Cold out
      const cnt_t2 = getAverage(['cnt-tc2-1', 'cnt-tc2-2', 'cnt-tc2-3']); // Cold in
      const cnt_Vh = parseFloat(document.getElementById('cnt-flow-h').value) || 2.0;
      const cnt_Vc = parseFloat(document.getElementById('cnt-flow-c').value) || 3.0;

      document.getElementById('lab-cnt-avg-T1').innerText = cnt_T1.toFixed(2) + ' ℃';
      document.getElementById('lab-cnt-avg-T2').innerText = cnt_T2.toFixed(2) + ' ℃';
      document.getElementById('lab-cnt-avg-t1').innerText = cnt_t1.toFixed(2) + ' ℃';
      document.getElementById('lab-cnt-avg-t2').innerText = cnt_t2.toFixed(2) + ' ℃';

      // 2. Parallel Flow Averages
      const par_T1 = getAverage(['par-t1-1', 'par-t1-2', 'par-t1-3']);
      const par_T2 = getAverage(['par-t2-1', 'par-t2-2', 'par-t2-3']);
      const par_t1 = getAverage(['par-tc1-1', 'par-tc1-2', 'par-tc1-3']); // Cold out
      const par_t2 = getAverage(['par-tc2-1', 'par-tc2-2', 'par-tc2-3']); // Cold in
      const par_Vh = parseFloat(document.getElementById('par-flow-h').value) || 2.0;
      const par_Vc = parseFloat(document.getElementById('par-flow-c').value) || 3.0;

      document.getElementById('lab-par-avg-T1').innerText = par_T1.toFixed(2) + ' ℃';
      document.getElementById('lab-par-avg-T2').innerText = par_T2.toFixed(2) + ' ℃';
      document.getElementById('lab-par-avg-t1').innerText = par_t1.toFixed(2) + ' ℃';
      document.getElementById('lab-par-avg-t2').innerText = par_t2.toFixed(2) + ' ℃';

      // Heat transfer rates: mass flow rates (kg/s)
      // V (L/min) * 1 kg/L / 60
      const cnt_mh = (cnt_Vh * CONSTANTS.rho) / 60.0;
      const cnt_mc = (cnt_Vc * CONSTANTS.rho) / 60.0;
      const cnt_qh = cnt_mh * CONSTANTS.Cp * (cnt_T1 - cnt_T2);
      const cnt_qc = cnt_mc * CONSTANTS.Cp * (cnt_t1 - cnt_t2);
      const cnt_qavg = (cnt_qh + cnt_qc) / 2.0;

      const par_mh = (par_Vh * CONSTANTS.rho) / 60.0;
      const par_mc = (par_Vc * CONSTANTS.rho) / 60.0;
      const par_qh = par_mh * CONSTANTS.Cp * (par_T1 - par_T2);
      const par_qc = par_mc * CONSTANTS.Cp * (par_t1 - par_t2);
      const par_qavg = (par_qh + par_qc) / 2.0;

      // LMTD calculation
      // Counterflow: dt1 = T1 - t1, dt2 = T2 - t2
      const cnt_dt1 = cnt_T1 - cnt_t1;
      const cnt_dt2 = cnt_T2 - cnt_t2;
      const cnt_lmtd = (Math.abs(cnt_dt1 - cnt_dt2) < 1e-4) ? cnt_dt1 : (cnt_dt1 - cnt_dt2) / Math.log(cnt_dt1 / cnt_dt2);

      // Parallel: dt1 = T1 - t2, dt2 = T2 - t1
      const par_dt1 = par_T1 - par_t2;
      const par_dt2 = par_T2 - par_t1;
      const par_lmtd = (Math.abs(par_dt1 - par_dt2) < 1e-4) ? par_dt1 : (par_dt1 - par_dt2) / Math.log(par_dt1 / par_dt2);

      // Overall Heat Transfer Coefficients
      const cnt_Ui = cnt_qavg / (CONSTANTS.Ai * cnt_lmtd);
      const cnt_Uo = cnt_qavg / (CONSTANTS.Ao * cnt_lmtd);
      const cnt_Um = cnt_qavg / (CONSTANTS.Am * cnt_lmtd);

      const par_Ui = par_qavg / (CONSTANTS.Ai * par_lmtd);
      const par_Uo = par_qavg / (CONSTANTS.Ao * par_lmtd);
      const par_Um = par_qavg / (CONSTANTS.Am * par_lmtd);

      // Update Summary Table
      document.getElementById('res-cnt-qh').innerText = Math.round(cnt_qh).toLocaleString() + ' W';
      document.getElementById('res-par-qh').innerText = Math.round(par_qh).toLocaleString() + ' W';
      document.getElementById('res-cnt-qc').innerText = Math.round(cnt_qc).toLocaleString() + ' W';
      document.getElementById('res-par-qc').innerText = Math.round(par_qc).toLocaleString() + ' W';
      document.getElementById('res-cnt-qavg').innerText = Math.round(cnt_qavg).toLocaleString() + ' W';
      document.getElementById('res-par-qavg').innerText = Math.round(par_qavg).toLocaleString() + ' W';

      const diff = ((cnt_qavg - par_qavg) / par_qavg) * 100;
      document.getElementById('res-q-ratio').innerText = `향류가 +${diff.toFixed(1)}% 우수`;

      document.getElementById('res-cnt-lmtd').innerText = cnt_lmtd.toFixed(2) + ' ℃';
      document.getElementById('res-par-lmtd').innerText = par_lmtd.toFixed(2) + ' ℃';

      document.getElementById('res-cnt-Ui').innerText = Math.round(cnt_Ui).toLocaleString() + ' W/m²·℃';
      document.getElementById('res-par-Ui').innerText = Math.round(par_Ui).toLocaleString() + ' W/m²·℃';
      document.getElementById('res-cnt-Uo').innerText = Math.round(cnt_Uo).toLocaleString() + ' W/m²·℃';
      document.getElementById('res-par-Uo').innerText = Math.round(par_Uo).toLocaleString() + ' W/m²·℃';
      document.getElementById('res-cnt-Um').innerText = Math.round(cnt_Um).toLocaleString() + ' W/m²·℃';
      document.getElementById('res-par-Um').innerText = Math.round(par_Um).toLocaleString() + ' W/m²·℃';

      // Render Lab Curves
      renderLabCharts(cnt_T1, cnt_T2, cnt_t1, cnt_t2, par_T1, par_T2, par_t1, par_t2);
    }

    function renderLabCharts(cnt_T1, cnt_T2, cnt_t1, cnt_t2, par_T1, par_T2, par_t1, par_t2) {
      // Counter chart: Hot (left to right: T1 -> T2), Cold (left to right: t1 -> t2)
      // At x=0 (Hot inlet): T=T1, Cold fluid is at outlet: T=t1
      // At x=L (Hot outlet): T=T2, Cold fluid is at inlet: T=t2
      if (!window.Chart) return;
      const canvasCnt = document.getElementById('chart-lab-counter');
      const canvasPar = document.getElementById('chart-lab-parallel');
      if (!canvasCnt || !canvasPar) return;

      const xLabels = ['입구 (A지점)', '전열 1/4', '중간 (1/2)', '전열 3/4', '출구 (B지점)'];
      
      // Interpolate curves smoothly
      const cntHot = [cnt_T1, cnt_T1*0.75 + cnt_T2*0.25, cnt_T1*0.5 + cnt_T2*0.5, cnt_T1*0.25 + cnt_T2*0.75, cnt_T2].map(v => v.toFixed(1));
      const cntCold = [cnt_t1, cnt_t1*0.75 + cnt_t2*0.25, cnt_t1*0.5 + cnt_t2*0.5, cnt_t1*0.25 + cnt_t2*0.75, cnt_t2].map(v => v.toFixed(1));

      const ctxCnt = canvasCnt.getContext('2d');
      if (chartLabCounter) {
        chartLabCounter.data.datasets[0].data = cntHot;
        chartLabCounter.data.datasets[1].data = cntCold;
        chartLabCounter.update('none');
      } else {
        chartLabCounter = new Chart(ctxCnt, {
          type: 'line',
          data: {
            labels: xLabels,
            datasets: [
              {
                label: '온수 실측 (고온)',
                data: cntHot,
                borderColor: '#ea580c', // Orange
                backgroundColor: 'rgba(234, 88, 12, 0.1)',
                borderWidth: 3,
                tension: 0.35,
                fill: false
              },
              {
                label: '냉수 실측 (저온)',
                data: cntCold,
                borderColor: '#2563eb', // Blue
                backgroundColor: 'rgba(37, 99, 235, 0.1)',
                borderWidth: 3,
                tension: 0.35,
                fill: false
              }
            ]
          },
          options: getChartOptions('위치 (A지점 → B지점)', '온도 (℃)')
        });
      }

      // Parallel chart: Hot (T1 -> T2), Cold (t2 -> t1)
      const parHot = [par_T1, par_T1*0.65 + par_T2*0.35, par_T1*0.4 + par_T2*0.6, par_T1*0.18 + par_T2*0.82, par_T2].map(v => v.toFixed(1));
      const parCold = [par_t2, par_t2*0.65 + par_t1*0.35, par_t2*0.4 + par_t1*0.6, par_t2*0.18 + par_t1*0.82, par_t1].map(v => v.toFixed(1));

      const ctxPar = canvasPar.getContext('2d');
      if (chartLabParallel) {
        chartLabParallel.data.datasets[0].data = parHot;
        chartLabParallel.data.datasets[1].data = parCold;
        chartLabParallel.update('none');
      } else {
        chartLabParallel = new Chart(ctxPar, {
          type: 'line',
          data: {
            labels: xLabels,
            datasets: [
              {
                label: '온수 실측 (고온)',
                data: parHot,
                borderColor: '#ea580c', // Orange
                backgroundColor: 'rgba(234, 88, 12, 0.1)',
                borderWidth: 3,
                tension: 0.35,
                fill: false
              },
              {
                label: '냉수 실측 (저온)',
                data: parCold,
                borderColor: '#2563eb', // Blue
                backgroundColor: 'rgba(37, 99, 235, 0.1)',
                borderWidth: 3,
                tension: 0.35,
                fill: false
              }
            ]
          },
          options: getChartOptions('위치 (A지점 → B지점)', '온도 (℃)')
        });
      }
    }

    function loadSampleLabData() {
      // Counterflow sample
      document.getElementById('cnt-flow-h').value = "2.0";
      document.getElementById('cnt-flow-c').value = "3.0";
      document.getElementById('cnt-t1-1').value = "58.2";
      document.getElementById('cnt-t2-1').value = "44.1";
      document.getElementById('cnt-tc1-1').value = "24.5";
      document.getElementById('cnt-tc2-1').value = "15.1";

      document.getElementById('cnt-t1-2').value = "58.0";
      document.getElementById('cnt-t2-2').value = "43.8";
      document.getElementById('cnt-tc1-2').value = "24.6";
      document.getElementById('cnt-tc2-2').value = "15.2";

      document.getElementById('cnt-t1-3').value = "57.8";
      document.getElementById('cnt-t2-3').value = "43.9";
      document.getElementById('cnt-tc1-3').value = "24.4";
      document.getElementById('cnt-tc2-3').value = "15.0";

      // Parallel flow sample
      document.getElementById('par-flow-h').value = "2.0";
      document.getElementById('par-flow-c').value = "3.0";
      document.getElementById('par-t1-1').value = "58.1";
      document.getElementById('par-t2-1').value = "46.5";
      document.getElementById('par-tc1-1').value = "22.7";
      document.getElementById('par-tc2-1').value = "15.1";

      document.getElementById('par-t1-2').value = "57.9";
      document.getElementById('par-t2-2').value = "46.4";
      document.getElementById('par-tc1-2').value = "22.8";
      document.getElementById('par-tc2-2').value = "15.0";

      document.getElementById('par-t1-3').value = "58.0";
      document.getElementById('par-t2-3').value = "46.6";
      document.getElementById('par-tc1-3').value = "22.6";
      document.getElementById('par-tc2-3').value = "15.2";

      calcLabResults();
    }

    function resetLabData() {
      ['cnt-t1-1', 'cnt-t2-1', 'cnt-tc1-1', 'cnt-tc2-1',
       'cnt-t1-2', 'cnt-t2-2', 'cnt-tc1-2', 'cnt-tc2-2',
       'cnt-t1-3', 'cnt-t2-3', 'cnt-tc1-3', 'cnt-tc2-3',
       'par-t1-1', 'par-t2-1', 'par-tc1-1', 'par-tc2-1',
       'par-t1-2', 'par-t2-2', 'par-tc1-2', 'par-tc2-2',
       'par-t1-3', 'par-t2-3', 'par-tc1-3', 'par-tc2-3'].forEach(id => {
        document.getElementById(id).value = "0.0";
      });
      calcLabResults();
    }

    function downloadCSVReport() {
      const rows = [
        ["구분", "모드", "회차", "고온유량(L/min)", "저온유량(L/min)", "T1(℃)", "T2(℃)", "t1(℃)", "t2(℃)"],
        ["실험데이터", "향류", "1회", document.getElementById('cnt-flow-h').value, document.getElementById('cnt-flow-c').value, document.getElementById('cnt-t1-1').value, document.getElementById('cnt-t2-1').value, document.getElementById('cnt-tc1-1').value, document.getElementById('cnt-tc2-1').value],
        ["실험데이터", "향류", "2회", document.getElementById('cnt-flow-h').value, document.getElementById('cnt-flow-c').value, document.getElementById('cnt-t1-2').value, document.getElementById('cnt-t2-2').value, document.getElementById('cnt-tc1-2').value, document.getElementById('cnt-tc2-2').value],
        ["실험데이터", "향류", "3회", document.getElementById('cnt-flow-h').value, document.getElementById('cnt-flow-c').value, document.getElementById('cnt-t1-3').value, document.getElementById('cnt-t2-3').value, document.getElementById('cnt-tc1-3').value, document.getElementById('cnt-tc2-3').value],
        ["실험데이터", "병류", "1회", document.getElementById('par-flow-h').value, document.getElementById('par-flow-c').value, document.getElementById('par-t1-1').value, document.getElementById('par-t2-1').value, document.getElementById('par-tc1-1').value, document.getElementById('par-tc2-1').value],
        ["실험데이터", "병류", "2회", document.getElementById('par-flow-h').value, document.getElementById('par-flow-c').value, document.getElementById('par-t1-2').value, document.getElementById('par-t2-2').value, document.getElementById('par-tc1-2').value, document.getElementById('par-tc2-2').value],
        ["실험데이터", "병류", "3회", document.getElementById('par-flow-h').value, document.getElementById('par-flow-c').value, document.getElementById('par-t1-3').value, document.getElementById('par-t2-3').value, document.getElementById('par-tc1-3').value, document.getElementById('par-tc2-3').value]
      ];

      let csvContent = "\uFEFF"; // UTF-8 BOM
      rows.forEach(r => {
        csvContent += r.join(",") + "\r\n";
      });

      const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.setAttribute("href", url);
      link.setAttribute("download", "이중관식열교환기_실험결과기록.csv");
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }

    // ==========================================
    // QUIZ CHECKER
    // ==========================================
    function checkAnswer(qName, correctAnswer, explanation) {
      const selected = document.querySelector(`input[name="${qName}"]:checked`);
      const feedbackEl = document.getElementById(`feedback-${qName}`);
      feedbackEl.classList.remove('hidden');

      if (!selected) {
        feedbackEl.innerHTML = `<span class="text-amber-600 dark:text-amber-400">보기를 선택한 후 정답 확인을 눌러주세요.</span>`;
        return;
      }

      if (selected.value === correctAnswer) {
        feedbackEl.innerHTML = `
          <div class="p-3 mt-2 rounded-xl bg-emerald-50 dark:bg-emerald-950/40 border border-emerald-300 dark:border-emerald-800 text-emerald-800 dark:text-emerald-300">
            <div class="font-bold flex items-center gap-1.5 mb-1"><i data-lucide="check" class="w-4 h-4"></i> 정답입니다!</div>
            <p class="text-[11px] leading-relaxed font-normal">${explanation}</p>
          </div>
        `;
      } else {
        feedbackEl.innerHTML = `
          <div class="p-3 mt-2 rounded-xl bg-rose-50 dark:bg-rose-950/40 border border-rose-300 dark:border-rose-800 text-rose-800 dark:text-rose-300">
            <div class="font-bold flex items-center gap-1.5 mb-1"><i data-lucide="x" class="w-4 h-4"></i> 다시 생각해보세요!</div>
            <p class="text-[11px] leading-relaxed font-normal">힌트: 열교환기의 국소 온도 구동력($\\Delta T$) 및 엔트로피/열역학 제2법칙을 고려해보세요.</p>
          </div>
        `;
      }
      lucide.createIcons();
    }

    // ==========================================
    // AI TUTOR & CHEMICAL ENGINEERING KNOWLEDGE ENGINE
    // ==========================================
    function toggleFloatingChat() {
      const win = document.getElementById('floating-chat-window');
      win.classList.toggle('hidden');
      if (!win.classList.contains('hidden')) {
        document.getElementById('float-chat-input').focus();
        lucide.createIcons();
      }
    }

    function openApiConfigModal() {
      document.getElementById('api-modal').classList.remove('hidden');
      const savedKey = localStorage.getItem('gemini_api_key') || '';
      document.getElementById('input-gemini-key').value = savedKey;
      lucide.createIcons();
    }

    function closeApiConfigModal() {
      document.getElementById('api-modal').classList.add('hidden');
    }

    function saveApiKey() {
      const key = document.getElementById('input-gemini-key').value.trim();
      if (key) {
        localStorage.setItem('gemini_api_key', key);
        alert('Gemini API 키가 저장되었습니다. 실시간 Gemini 모델로 질의응답이 연동됩니다.');
      }
      closeApiConfigModal();
      updateApiBadges();
    }

    function clearApiKey() {
      localStorage.removeItem('gemini_api_key');
      document.getElementById('input-gemini-key').value = '';
      alert('API 키가 삭제되었습니다. 내장 화공 전문 지식엔진으로 전환됩니다.');
      closeApiConfigModal();
      updateApiBadges();
    }

    function updateApiBadges() {
      const hasKey = !!localStorage.getItem('gemini_api_key');
      const dot = document.getElementById('api-status-badge');
      if (dot) {
        dot.className = hasKey ? 'w-2 h-2 rounded-full bg-emerald-500' : 'w-2 h-2 rounded-full bg-slate-400';
      }
      const floatStatus = document.getElementById('float-api-status');
      if (floatStatus) {
        floatStatus.innerText = hasKey ? 'Gemini 1.5/2.5 Live 연동 중' : '화공 전문 지식엔진 동작 중';
      }
    }

    // Built-in Knowledge Base & Rule-based NLP Engine
    const KNOWLEDGE_BASE = [
      {
        keywords: ['향류', '병류', '효율', '비교', '차이점', '왜 향류'],
        title: '향류(Counterflow)가 병류(Parallel Flow)보다 효율적인 이유',
        response: `**향류(Counterflow)가 병류보다 항상 높은 열효율을 나타내는 핵심 원리:**

1. **전 구간 균일한 전열 구동력($\\Delta T$):**
   - **병류**는 입구에서 온도차가 매우 크지만, 출구로 갈수록 온수와 냉수의 온도가 가까워지면서 열전달 구동력이 급격히 $0$에 수렴합니다.
   - 반면 **향류**는 고온 입구와 냉수 출구, 온수 출구와 냉수 입구가 각각 맞물려 장치 전 구간에 걸쳐 일정한 $\\Delta T$ 구동력을 유지합니다.
2. **높은 대수평균온도차($\\Delta t_m$):**
   - 열전달 기본식 $q = U A \\Delta t_m$에서 동일한 열전달 면적($A$)과 총괄계수($U$)일 때, 향류의 LMTD가 병류의 LMTD보다 항상 크므로 더 많은 열량($q$)을 전달할 수 있습니다.
3. **온도 교차(Temperature Cross) 가능:**
   - 향류에서는 냉각수 출구 온도가 온수 출구 온도보다 높아지는 것이 가능하여 고온 유체의 열을 훨씬 깊이 회수할 수 있습니다.`
      },
      {
        keywords: ['lmtd', '대수평균', '산술평균', '로그', '유도', '왜 lmtd'],
        title: 'LMTD(대수평균온도차)의 수학적 유도 및 산술평균을 쓰지 않는 이유',
        response: `**대수평균온도차(LMTD, $\\Delta t_m$)를 반드시 사용해야 하는 이유:**

- 유체가 관을 따라 흐를 때 방출/흡수되는 미소 열량은 $dq = -m_h C_{ph} dT_h = m_c C_{pc} dt_c = U \\Delta T dA$ 입니다.
- 온도차의 변화율을 미분방정식으로 전개하면:
  $$\\frac{d(\\Delta T)}{dA} = -U \\left(\\frac{1}{m_h C_{ph}} + \\frac{1}{m_c C_{pc}}\\right) \\Delta T$$
- 이는 온도차가 거리에 따라 선형이 아니라 **지수함수($e^{-kx}$)** 형태로 감소함을 의미합니다.
- 이를 전열 면적 $A$에 대해 적분하면:
  $$q = U A \\left( \\frac{\\Delta t_1 - \\Delta t_2}{\\ln(\\Delta t_1 / \\Delta t_2)} \\right)$$
- 따라서 유효 구동력은 산술평균이 아닌 **대수평균(LMTD)**이 됩니다. 산술평균을 사용하면 LMTD보다 항상 과대평가되어 열교환기를 실제보다 작게 설계하게 되는 심각한 문제가 발생합니다.`
      },
      {
        keywords: ['온도 교차', '교차', 'temperature cross', '역전'],
        title: '온도 교차(Temperature Cross)의 물리적 의미와 조건',
        response: `**온도 교차(Temperature Cross, $t_{c,out} > T_{h,out}$)의 의미:**

- **정의:** 차가운 냉각수의 출구 온도($t_1$)가 뜨거운 온수의 출구 온도($T_2$)보다 더 높아지는 현상입니다.
- **가능 여부:**
  - **향류:** 가능함! 냉각수가 최종적으로 빠져나가는 곳은 온수의 가장 뜨거운 입구($T_1$) 지점과 열교환하므로, 냉각수 온도가 $T_2$를 넘어 $T_1$ 근처까지 상승할 수 있습니다.
  - **병류:** 불가능함! 두 유체가 같은 출구로 나가기 때문에 $t_1 = T_2$가 되는 순간 온도차가 $0$이 되어 열전달이 멈추며, 이를 넘어서면 열역학 제2법칙(차가운 곳에서 뜨거운 곳으로 비자발적 열이동)에 위배됩니다.
- 시뮬레이터 탭에서 '온도 교차 프리셋'을 누르시면 실시간으로 이 현상을 관찰하실 수 있습니다!`
      },
      {
        keywords: ['오차', '원인', '불일치', 'qh', 'qc', '열손실', '열수지'],
        title: '실험 시 온수 방출열량(qh)과 냉수 흡수열량(qc)의 오차 원인',
        response: `**화공실험2 이중관식 열교환기 실험 오차 원인 핵심 정리:**

1. **대기로의 외부 방열 손실 (가장 주된 원인):**
   - 외관(40A 동관)이 외기에 노출되어 있어 온수의 열 중 일부가 냉각수로 가지 않고 실험실 공기로 손실됩니다. 따라서 일반적으로 $q_h > q_c$가 됩니다.
2. **비정상 상태(Unsteady State) 계측:**
   - 밸브를 조작하여 유량을 바꾼 후 열적 평형에 도달하기까지 약 3~5분이 소요됩니다. 충분히 기다리지 않고 기록하면 큰 오차가 발생합니다.
3. **온도계 감온부(써모커플) 깊이 및 위치:**
   - 관내 유체는 유속 분포와 온도 구배(벽면 근처와 중심부 온도 차이)가 존재하므로 센서 접촉 깊이에 따라 편차가 생깁니다.
4. **보일러 온도 헌팅 & 유량계 눈금 판독 오차:**
   - 펌프 진동 및 눈금 시차(Parallax error)로 유량 측정에 오차가 개입됩니다.`
      },
      {
        keywords: ['면적', 'ai', 'ao', 'am', '표면적', '내경', '외경'],
        title: '총괄전열계수(U) 계산 시 어떤 면적(Ai, Ao, Am)을 써야 하나요?',
        response: `**전열면적($A_i, A_o, A_m$)과 총괄전열계수의 관계:**

- 열전달률은 일정하므로 $q = U_i A_i \\Delta t_m = U_o A_o \\Delta t_m = U_m A_m \\Delta t_m$ 관계가 성립합니다.
- **내표면적 기준 ($A_i = \\pi D_i L$):** 동관 내경(14.26 mm) 기준, $A_i \\approx 0.0672\\text{ m}^2$
- **외표면적 기준 ($A_o = \\pi D_o L$):** 동관 외경(15.88 mm) 기준, $A_o \\approx 0.0748\\text{ m}^2$
- **대수평균 표면적 ($A_m = \\frac{A_o - A_i}{\\ln(A_o / A_i)}$):** 관벽 두께(0.8 mm)를 고려한 엄밀한 평균 면적으로 $A_m \\approx 0.0709\\text{ m}^2$
- **권장 기준:** 일반적으로 관벽 두께가 얇은 박벽관(Thin-walled tube)에서는 외표면적($A_o$) 또는 대수평균면적($A_m$)을 주로 사용하며, 보고서에는 $U_i$와 $U_o$, $U_m$의 차이를 관벽 곡률 효과로 비교 설명하면 높은 점수를 받습니다.`
      },
      {
        keywords: ['고찰', '보고서', '작성법', '토의', 'discussion'],
        title: '화공실험2 열교환기 보고서 고찰(Discussion) 작성 가이드',
        response: `**보고서 고찰(Discussion) 작성 시 꼭 포함해야 할 A+ 체크리스트:**

1. **병류 vs 향류 실험 결과 정량 비교:**
   - 동일 유량 조건에서 산출된 $q_{avg}$, $\\Delta t_m$, $U$ 값을 표로 대조하고, 향류가 병류에 비해 몇 % 더 열전달 능력이 우수한지 수치로 기술.
2. **열수지 효율 (Heat Balance) 분석:**
   - 온수의 손실열($q_h$)과 냉각수의 흡수열($q_c$)의 비 $\\frac{q_c}{q_h} \\times 100(\\%)$를 구하고, 100%가 되지 않는 이유(대기 방열 손실 등)를 기술.
3. **유량 변화에 따른 $U$값 경향성 고찰:**
   - 유량이 증가할 때 유속 증가 $\\to$ 레이놀즈수($Re$) 증가 $\\to$ 난류 강화 $\\to$ 경계층 두께 감소 $\\to$ 대류계수 $h$ 상승 $\\to$ 총괄계수 $U$ 상승 논리 전개.
4. **Dittus-Boelter 식 이론치와의 오차율 비교:**
   - 문헌상의 이론적 대류 열전달계수와 본 실험의 실측 $U$값을 비교하고, 관내 스케일(Fouling factor) 및 센서 오차의 영향 언급.`
      },
      {
        keywords: ['유량', '속도', '증가', '레이놀즈', 'reynolds'],
        title: '유량 증가 시 열전달 특성 변화',
        response: `**유량이 증가할 때 일어나는 전열 특성:**

1. **총괄전열계수($U$) 증가:**
   - 유량이 늘어나면 유속 $v$가 빨라져 난류(Turbulent flow)가 발달합니다.
   - 열전달 저항의 대부분을 차지하는 관벽 접촉면의 '온도 경계층(Thermal boundary layer)'이 얇아져 대류 열전달계수 $h$가 급격히 상승합니다.
2. **출구 온도차의 변화:**
   - 유량이 많아지면 유체의 열용량 유량($C = m C_p$)이 커져 관을 통과하는 동안 유체의 온도 변화폭($\\Delta T$)은 오히려 줄어들지만, 총 열전달량($q = m C_p \\Delta T$)은 증가합니다.`
      }
    ];

    async function callGeminiAPI(apiKey, query) {
      const url = `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=${apiKey}`;
      const systemPrompt = `당신은 화학공학실험2 이중관식 열교환기(Double Pipe Heat Exchanger) 전문 AI 교수/조교입니다.
실험 규격: 내관 동관 15A (내경 14.26mm, 외경 15.88mm), 외관 40A, 전열길이 1500mm, 항온수조 58℃, 물 비열 4184 J/kg·K (0.999 cal/g·℃).
학생의 질문에 대해 학술적으로 정확하고 친절하며 논리정연하게 한국어로 답변해주세요. 수식은 KaTeX 기호($...$)를 사용하세요.`;

      const body = {
        contents: [
          {
            role: "user",
            parts: [{ text: `${systemPrompt}\n\n학생 질문: ${query}` }]
          }
        ]
      };

      const resp = await fetch(url, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body)
      });

      if (!resp.ok) {
        throw new Error(`API 오류: ${resp.status}`);
      }
      const data = await resp.json();
      return data.candidates[0].content.parts[0].text;
    }

    function searchKnowledgeBase(query) {
      const q = query.toLowerCase().replace(/[\s\?\,\.\!]/g, '');
      let bestMatch = null;
      let maxScore = 0;

      for (const item of KNOWLEDGE_BASE) {
        let score = 0;
        for (const kw of item.keywords) {
          const cleanKw = kw.toLowerCase().replace(/\s/g, '');
          if (q.includes(cleanKw)) {
            score += 3;
          }
        }
        if (score > maxScore) {
          maxScore = score;
          bestMatch = item;
        }
      }

      if (bestMatch && maxScore >= 3) {
        return bestMatch.response;
      }

      // Heuristic fallback response
      return `**화공실험 AI 튜터 해설:**

질문하신 내용: **"${query}"** 에 대한 핵심 안내입니다.

이중관식 열교환기(Double Pipe Heat Exchanger) 실험의 핵심 원리는 **내관의 고온 유체(온수)와 환상공간(Annulus)의 저온 유체(냉각수)** 사이의 전열 과정입니다.
- **열전달량 계산:** $q = m_h C_{ph} (T_1 - T_2) = m_c C_{pc} (t_1 - t_2)$
- **대수평균온도차:** $\\Delta t_m = \\frac{\\Delta t_1 - \\Delta t_2}{\\ln(\\Delta t_1 / \\Delta t_2)}$
- **총괄전열계수:** $U = \\frac{q}{A \\Delta t_m}$

만약 더 세부적인 자유 대화가 필요하시다면 우측 상단의 **[API 연동 설정]** 버튼을 통해 Google Gemini API Key를 등록하시면 실시간 LLM 답변을 무제한으로 받으실 수 있습니다!`;
    }

    async function handleUserQuery(query, outputContainerId) {
      const container = document.getElementById(outputContainerId);
      if (!container) return;

      // Append user message
      const userBubble = document.createElement('div');
      userBubble.className = "flex items-start justify-end gap-2";
      userBubble.innerHTML = `
        <div class="bg-purple-600 text-white p-3 rounded-2xl rounded-tr-sm max-w-xl shadow-sm text-xs leading-relaxed">
          ${escapeHTML(query)}
        </div>
        <div class="w-7 h-7 rounded-full bg-slate-200 dark:bg-slate-700 text-slate-700 dark:text-slate-300 flex items-center justify-center flex-shrink-0 font-bold text-[10px]">
          학생
        </div>
      `;
      container.appendChild(userBubble);
      container.scrollTop = container.scrollHeight;

      // Loading bubble
      const loadingBubble = document.createElement('div');
      loadingBubble.className = "flex items-start gap-2 ai-loading-indicator";
      loadingBubble.innerHTML = `
        <div class="w-7 h-7 rounded-full bg-gradient-to-tr from-purple-600 to-indigo-600 text-white flex items-center justify-center flex-shrink-0 text-[10px]">
          <i data-lucide="bot" class="w-3.5 h-3.5 animate-spin"></i>
        </div>
        <div class="bg-white dark:bg-slate-800 p-3 rounded-2xl rounded-tl-sm border border-slate-200 dark:border-slate-700 text-slate-500 text-xs italic">
          화학공학 지식엔진이 답변을 정리하고 있습니다...
        </div>
      `;
      container.appendChild(loadingBubble);
      container.scrollTop = container.scrollHeight;
      lucide.createIcons();

      let answer = "";
      const apiKey = localStorage.getItem('gemini_api_key');
      try {
        if (apiKey) {
          answer = await callGeminiAPI(apiKey, query);
        } else {
          // Simulate slight thoughtful delay for realism
          await new Promise(r => setTimeout(r, 350));
          answer = searchKnowledgeBase(query);
        }
      } catch (err) {
        console.error(err);
        answer = searchKnowledgeBase(query) + `\n\n*(참고: Gemini API 호출 오류로 내장 화공 지식엔진으로 자동 대체되었습니다.)*`;
      }

      // Remove loading indicator
      loadingBubble.remove();

      // Append AI Response Bubble
      const aiBubble = document.createElement('div');
      aiBubble.className = "flex items-start gap-2";
      aiBubble.innerHTML = `
        <div class="w-7 h-7 rounded-full bg-gradient-to-tr from-purple-600 to-indigo-600 text-white flex items-center justify-center flex-shrink-0 text-[10px] shadow">
          <i data-lucide="bot" class="w-3.5 h-3.5"></i>
        </div>
        <div class="bg-white dark:bg-slate-800 p-4 rounded-2xl rounded-tl-sm border border-slate-200 dark:border-slate-700 shadow-sm max-w-2xl text-slate-800 dark:text-slate-200 text-xs leading-relaxed space-y-2">
          ${formatMarkdownToHTML(answer)}
        </div>
      `;
      container.appendChild(aiBubble);
      container.scrollTop = container.scrollHeight;

      // Render KaTeX in the new bubble
      if (window.renderMathInElement) {
        renderMathInElement(aiBubble, {
          delimiters: [
            {left: '$$', right: '$$', display: true},
            {left: '$', right: '$', display: false}
          ],
          throwOnError: false
        });
      }
      lucide.createIcons();
    }

    function escapeHTML(str) {
      return str.replace(/[&<>'"]/g, tag => ({
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        "'": '&#39;',
        '"': '&quot;'
      }[tag] || tag));
    }

    function formatMarkdownToHTML(text) {
      // Basic markdown conversion (bold, lists, code, paragraphs)
      let html = text
        .replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>')
        .replace(/\*(.*?)\*/g, '<em>$1</em>')
        .replace(/`([^`]+)`/g, '<code class="px-1.5 py-0.5 rounded bg-slate-100 dark:bg-slate-700 text-purple-600 dark:text-purple-400 font-mono text-[11px]">$1</code>');
      
      const paragraphs = html.split('\n\n');
      return paragraphs.map(p => {
        if (p.trim().startsWith('- ') || p.trim().startsWith('1. ')) {
          const lines = p.split('\n');
          let listHtml = '<ul class="list-disc list-inside space-y-1 my-1">';
          lines.forEach(l => {
            const clean = l.replace(/^[-*]\s+|\d+\.\s+/, '');
            listHtml += `<li>${clean}</li>`;
          });
          listHtml += '</ul>';
          return listHtml;
        }
        return `<p>${p.replace(/\n/g, '<br>')}</p>`;
      }).join('');
    }

    // Handlers for Tab and Floating chat
    function submitAIQuery() {
      const input = document.getElementById('ai-user-query');
      const val = input.value.trim();
      if (!val) return;
      input.value = '';
      handleUserQuery(val, 'ai-chat-output');
    }

    function submitFloatQuery() {
      const input = document.getElementById('float-chat-input');
      const val = input.value.trim();
      if (!val) return;
      input.value = '';
      handleUserQuery(val, 'float-chat-messages');
    }

    function askAITutor(question) {
      // Also switch to QA tab or open floating chat if not already there
      const currentTab = document.querySelector('.tab-content:not(.hidden)').id;
      if (currentTab !== 'tab-qa') {
        const floatWin = document.getElementById('floating-chat-window');
        if (floatWin.classList.contains('hidden')) {
          toggleFloatingChat();
        }
        handleUserQuery(question, 'float-chat-messages');
      } else {
        handleUserQuery(question, 'ai-chat-output');
      }
    }

    function clearAIChat() {
      document.getElementById('ai-chat-output').innerHTML = `
        <div class="flex items-start gap-3">
          <div class="w-8 h-8 rounded-full bg-gradient-to-tr from-purple-600 to-indigo-600 text-white flex items-center justify-center flex-shrink-0 font-bold text-xs shadow">
            <i data-lucide="bot" class="w-4 h-4"></i>
          </div>
          <div class="bg-white dark:bg-slate-800 p-4 rounded-2xl rounded-tl-sm border border-slate-200 dark:border-slate-700 shadow-sm max-w-2xl text-slate-800 dark:text-slate-200 leading-relaxed space-y-2">
            <p class="font-bold text-purple-600 dark:text-purple-400">대화가 초기화되었습니다. 새로운 질문을 입력해 주세요! 🧪</p>
          </div>
        </div>
      `;
      lucide.createIcons();
    }

    // Initial check for API badges & initial computations on load
    window.addEventListener('DOMContentLoaded', () => {
      try { updateApiBadges(); } catch (e) { console.warn(e); }
      try { updateSimulation(); } catch (e) { console.warn(e); }
      try { calcLabResults(); } catch (e) { console.warn(e); }
      if (window.lucide && lucide.createIcons) {
        lucide.createIcons();
      }
    });
