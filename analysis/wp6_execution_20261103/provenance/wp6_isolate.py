from pathlib import Path
import shutil,json,hashlib
root=Path.cwd();src=root/'model/transport_v1/source';dst=root/'model/transport_wp6_v1/source'
assert not dst.exists();shutil.copytree(src,dst)
p=dst/'appl/simpleRichardsnc/transportdiagnostics.hh';s=p.read_text()
s=s.replace('std::vector<double> dofZ_;','std::vector<double> dofZ_, dofX_;\n    std::vector<std::size_t> dofLevel_;\n    std::size_t levelCount_ = 0;')
s=s.replace('std::ofstream budget_, surface_;','std::ofstream budget_, surface_, weighted_;')
s=s.replace('"stop_record.json", "transport_metadata.json"','"stop_record.json", "transport_metadata.json", "weighted_online.csv"')
s=s.replace('dofZ_.resize(gg.numDofs());','dofZ_.resize(gg.numDofs()); dofX_.resize(gg.numDofs()); dofLevel_.resize(gg.numDofs());')
s=s.replace('levels.push_back(z); dofZ_[gg.vertexMapper().index(v)] = z;','levels.push_back(z); dofZ_[gg.vertexMapper().index(v)] = z;\n            dofX_[gg.vertexMapper().index(v)] = v.geometry().center()[0];')
s=s.replace('for (double fraction : {0.9, 0.75, 0.5}) {','''levelCount_ = levels.size();
        for (std::size_t d=0; d<dofZ_.size(); ++d)
            dofLevel_[d] = std::lower_bound(levels.begin(), levels.end(), dofZ_[d])-levels.begin();
        const auto lockedCuts = getParam<std::vector<double>>("Transport.CutCoordinates");
        if (lockedCuts.size()!=7) throw std::runtime_error("WP6 requires seven locked cut coordinates");
        std::size_t cutIndex=0;
        for (double fraction : {0.9, 0.75, 0.5, 0.995, 0.99, 0.975, 0.95}) {''')
s=s.replace('const double target = bottom_+fraction*height_;','const double target = lockedCuts[cutIndex++];')
s=s.replace('regions_.push_back({"upper_"+std::to_string(fraction), fraction, best});','if (std::abs(best-target)>1e-12) throw std::runtime_error("Locked WP6 cut is not a conservative BOX face");\n            regions_.push_back({"upper_"+std::to_string(fraction), fraction, best});')
s=s.replace('budget_.open("transport_budget.csv"); surface_.open("surface_events.csv");','''budget_.open("transport_budget.csv"); surface_.open("surface_events.csv");
        weighted_.open("weighted_online.csv"); weighted_ << std::setprecision(17);
        weighted_ << "step,time_s,weighted_top_mean,weighted_top_sigma,weighted_fourier,weighted_projection,volume_perturbation_rms,harmonic2_top,harmonic3_top,projection_defined\\n";''')
s=s.replace('schema_version\\\":1,\\\"inventory','schema_version\\\":2,\\\"inventory')
s=s.replace('double topWeighted=0, topArea=0;','double topWeighted=0, topArea=0;\n        std::vector<double> topWeight(x.size(),0.0), volumeWeight(x.size(),0.0);')
s=s.replace('const double z=dofZ_[scv.dofIndex()];','volumeWeight[scv.dofIndex()]+=scv.volume()/width_;\n                const double z=dofZ_[scv.dofIndex()];')
s=s.replace('topArea+=face.area(); topWeighted+=face.area()*x[inside.dofIndex()][1];','topArea+=face.area(); topWeighted+=face.area()*x[inside.dofIndex()][1];\n                        topWeight[inside.dofIndex()]+=face.area();')
needle='        budget_.flush(); surface_.flush(); checkStream(budget_); checkStream(surface_);'
block='''        // Read-only weighted statistics; never enter the physical residual.
        const double mean=topWeighted/topArea;
        std::vector<double> rowSum(levelCount_,0.0), rowWeight(levelCount_,0.0);
        for (std::size_t d=0;d<x.size();++d) {
            rowSum[dofLevel_[d]]+=volumeWeight[d]*x[d][1];
            rowWeight[dofLevel_[d]]+=volumeWeight[d];
        }
        for (std::size_t l=0;l<levelCount_;++l) rowSum[l]/=rowWeight[l];
        double wvar=0, re=0, im=0, proj=0, den=0, rms=0, volume=0;
        std::array<double,2> hre={0,0}, him={0,0};
        const double k=2*std::acos(-1.0)*problem_.waveNumber()/width_;
        for (std::size_t d=0;d<x.size();++d) {
            const double phase=k*(dofX_[d]-gg_.bBoxMin()[0]);
            const double v=x[d][1]-mean;
            wvar+=topWeight[d]*v*v; re+=topWeight[d]*v*std::cos(phase); im-=topWeight[d]*v*std::sin(phase);
            for (std::size_t j=0;j<2;++j) {
                hre[j]+=topWeight[d]*v*std::cos((j+2)*phase);
                him[j]-=topWeight[d]*v*std::sin((j+2)*phase);
            }
            const double anomaly=x[d][1]-rowSum[dofLevel_[d]];
            rms+=volumeWeight[d]*anomaly*anomaly; volume+=volumeWeight[d];
            if (problem_.hasEigenprofile()) {
                const double basis=problem_.eigenAmplitudeAt(dofZ_[d])*std::cos(phase);
                proj+=volumeWeight[d]*anomaly*basis; den+=volumeWeight[d]*basis*basis;
            }
        }
        weighted_ << step_ << ',' << time << ',' << mean << ',' << std::sqrt(wvar/topArea)
                  << ',' << 2*std::hypot(re,im)/width_ << ',';
        if (den>0) weighted_ << proj/den; // Background projection is undefined, not zero.
        weighted_ << ',' << std::sqrt(rms/volume)
                  << ',' << 2*std::hypot(hre[0],him[0])/width_ << ',' << 2*std::hypot(hre[1],him[1])/width_
                  << ',' << int(den>0) << '\\n';
        weighted_.flush(); checkStream(weighted_);
'''+needle
assert needle in s;s=s.replace(needle,block);p.write_text(s,encoding='utf-8',newline='\n')
manifest={str(q.relative_to(dst)).replace('\\','/'):hashlib.sha256(q.read_bytes()).hexdigest() for q in dst.rglob('*') if q.is_file()}
(dst.parent/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for name in ['main.cc','model.hh']:assert (dst/'appl/simpleRichardsnc'/name).read_bytes()==(src/'appl/simpleRichardsnc'/name).read_bytes()
print('Isolated source prepared; physics main.cc/model.hh byte-identical.')
