var gulp = require('gulp');
var jscs = require('gulp-jscs');

gulp.task('lint', () => {
    return gulp.src('src/**/*.js')
        .pipe(jscs())
        .pipe(jscs.reporter());
});
